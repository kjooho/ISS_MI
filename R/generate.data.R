##################################################################
#              Functions to Generate Case-Cohort Data            #
##################################################################
# Simulation data: 4 low-cost covariates and 2 expensive biomarkers (one continuous, one binary)
simul.data <- function(n, max.t = 15, Smax.t = 0.98, 
                       alpha = 1, relz1 = exp(1.5), relz2 = exp(0.5), relz3 = exp(0.1),
                       relxc = exp(0.1), relxb = exp(0.2),
                       relz1.xc = NULL,
                       max.enter = 2, Smax.censor = 0.9, 
                       mu  = c(z0 = 0, z1 = 0),
                       sig = c(sd0 = 1, sd1 = 1),
                       rho = c(rho01 = 0.05), 
                       p2 = 0.5, p3 = 0.3,
                       r_obs = c(r02 = -0.05, r12 = 0.01),
                       side = c("lt","gt"),
                       xc_coef, xb_coef,
                       missing.frac, design = "case.cohort") {
  z0 <- rnorm(n, mean = mu["z0"], sd = sig["sd0"])
  z1 <- rnorm(n, mean = mu["z1"], sd = sig["sd1"])
  z2 <- rbinom(n, size = 1, prob = p2)
  z3 <- rbinom(n, size = 1, prob = p3)
  
  # Continuous expensive covariates
  make_xc <- function(coef, err) {
    with(as.list(coef), a1*z0 + a2*z1 + a3*z2 + a4*z3 + err)
  }
  
  make_xb <- function(coef) {
    lin_xb <- with(as.list(coef), b1*z0 + b2*z1 + b3*z2 + b4*z3)
    p_xb <- 1 / (1 + exp(-lin_xb))
    xb <- rbinom(n, 1, p_xb)
    return(xb)
  }
  xc <- make_xc(xc_coef, rnorm(n, 0, 1))
  xc <- (xc - mean(xc)) / sd(xc)
  xb <- make_xb(xb_coef)
  
  # Quadratic and interaction variables
  z1 <- (z1 - mean(z1)) / sd(z1)
  z1.xc <- z1 * xc
  
  # coefficients
  b_z0 <- log(-log(Smax.t) / (max.t^alpha))
  b_z1 <- log(relz1)
  b_z2 <- log(relz2)
  b_z3 <- log(relz3)
  b_xc <- log(relxc)
  b_xb <- log(relxb)
  
  eta <- b_z0 + b_z1 * z1 + b_z2 * z2 + b_z3 * z3 + b_xc*xc + b_xb * xb
  
  if (!is.null(relz1.xc))  eta <- eta + log(relz1.xc) * z1.xc
  lambda <- exp(eta / alpha)
  
  # survival data
  survtime = stats::rweibull(n, shape = alpha, scale = 1 / lambda)
  entrtime = stats::runif(n, 0, max.enter)
  deadtime = stats::rexp(n, rate = -log(Smax.censor)/max.t)
  censtime = max.t - entrtime
  eventime = pmin(survtime, censtime, deadtime)
  ind.fail = as.integer(survtime <= pmin(censtime, deadtime))
  
  if (is.null(relz1.xc)) {covariates <- data.frame(z0, z1, z2, z3, xc, xb)
  } else if (!is.null(relz1.xc)) {covariates <- data.frame(z0, z1, z2, z3, xc, xb, z1.xc)
  } 
  DATA = data.frame(cbind(covariates, eventime, ind.fail))
  N <- nrow(DATA)
  if (design == "stratified.case.cohort"){
    if(!require(sampling)){
      install.packages("sampling")
      library(sampling)
    }
    # Stratified case-cohort
    DATA$stratum <- interaction(z3, DATA$z2, drop = TRUE)
    tbl  <- table(DATA$stratum)
    levs <- names(tbl)
    prop <- as.numeric(tbl) / sum(tbl)
    l <- round(N * (1 - missing.frac))
    size_raw <- l * prop
    size <- floor(size_raw)
    k <- l - sum(size)
    if (k > 0) {
      add <- order(size_raw - size, decreasing = TRUE)[seq_len(k)]
      size[add] <- size[add] + 1L
    }
    DATA <- DATA[order(DATA[["stratum"]]), ] # Must arrange the data by their strata!
    init.sample <- as.integer(sampling::strata(DATA, stratanames = "stratum", 
                                               size = size, method = "srswor")$ID_unit) # MUST specify "sampling::"strata
  } else if (design == "case.cohort"){
    init.sample <- sample(1:N, round(N * (1 - missing.frac)))
  } else {
    stop("design must be 'stratified.case.cohort' or 'case.cohort'.")
  }
  leftover <- setdiff(1:N, init.sample)
  leftover.cases <- leftover[DATA$ind.fail[leftover] == 1]
  final.sample <- c(init.sample, leftover.cases)
  DATA$ind.ph2 <- ifelse((1:N) %in% final.sample, 1, 0)
  DATA$ind.ini <- ifelse((1:N) %in% init.sample, 1, 0)
  DATA$id <- 1:N
  DATA$weights <- 0
  return(DATA)
}

##################################################################
# generate.data()
#
# Inputs:
#   missing.frac : Proportion of units with missing expensive biomarkers.
#   design       : Sampling design indicator (e.g., case-cohort, supersample, etc.).
#   interaction  : Logical; if TRUE, include an interaction term (z1 × xc1).
#
# Details:
#   - Calls simul.data() with fixed baseline parameters (n=25,000, follow-up=15).
#   - Optionally includes an interaction effect through relz1.xc1.
#   - Converts xb to a factor and imposes phase-2 missingness on xc, xb, and
#     the interaction term.
#   - Computes baseline cumulative hazard using a null Cox model to support
#     SMC-FCS compatibility.
#
# Returns:
#   A list with:
#     dat     : Full simulated dataset (complete data).
#     dat.mi  : Dataset with imposed missingness for MI.
#     dat.smc : Dataset formatted for SMC-FCS imputation.
#
##################################################################
generate.data <- function(missing.frac, design, interaction){
  if (!is.logical(interaction) || length(interaction) != 1) {
    stop("'interaction' must be a single logical value (TRUE or FALSE).")
  }
  args_common <- list(n = 25000, max.t = 15, Smax.t = 0.998, 
    alpha = 1, relz1 = exp(1.5), relz2 = exp(0.5),
    relxc = exp(0.4), relxb = exp(0.3),
    max.enter = 2, Smax.censor = 0.9, 
    mu  = c(z0 = 0, z1 = 0),
    sig = c(sd0 = 1, sd1 = 1),
    rho = c(rho01 = 0.05), 
    p2 = 0.5,
    r_obs = c(r02 = -0.05, r12 = 0.01),
    side = c("lt","gt"),
    xc_coef = c(a1=0.20, a2=0.10, a3=0.10, a4=0.07),
    xb_coef = c(b1=0.15, b2=0.10, b3=0.07, b4=0.05),
    missing.frac = missing.frac, design = design)
  
  if (interaction) {
    args_common$relz1.xc1 <- exp(0.3)
    args_common$Smax.t <- 0.9989
  }
  
  dat <- do.call(simul.data, args_common)
  
  dat$xb <- factor(dat$xb)
  cols <- if (interaction) c("xc", "xb", "z1.xc") else c("xc", "xb")
  dat.smc <- dat
  dat.smc[, cols] <- lapply(dat.smc[, cols, drop = FALSE], function(v) {
    v[dat$ind.ph2 == 0] <- NA
    v
  })
  dat.mi <- dat.smc
  
  fits <- coxph(Surv(eventime, ind.fail) ~ 1, data = dat.mi)
  sf <- survfit(fits)
  na.df <- data.frame(time = sf$time, cumhaz = sf$cumhaz)
  na.df <- na.df[order(na.df$time), , drop = FALSE]
  
  # Find baseline hazard in case of missing values in `cumhaz` in `na.df`
  idx <- findInterval(dat.mi$eventime, na.df$time, left.open = FALSE)
  dat.mi$baseline_hazard <- ifelse(idx == 0, 0, na.df$cumhaz[idx])
  return(list(
    dat = dat,
    dat.mi = dat.mi,
    dat.smc = dat.smc
  ))
}


