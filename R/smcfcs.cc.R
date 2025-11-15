smcfcs.cc <- function (originaldata, smtype, smformula, predictorMatrix = NULL, method,
                       m = 5, numit = 10, rjlimit = 1000, noisy = FALSE, errorProneMatrix = NULL, in.subco, super) 
{
  smcfcs.core.cc(originaldata, smtype, smformula, predictorMatrix, method,
                 m, numit, rjlimit, noisy, errorProneMatrix = errorProneMatrix, in.subco, super)
}


smcfcs.core.cc <- function (originaldata, smtype, smformula, predictorMatrix = NULL, method,
                            m = 5, numit = 10, rjlimit = 1000, noisy = FALSE, errorProneMatrix = NULL, in.subco, super,
                            ...) 
{
  extraArgs <- list(...)
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata) != length(method)) 
    stop("Method argument must have the same length as the number of columns in the data frame.")
  n <- dim(originaldata)[1]
  
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)
  
  if ((smtype %in% c("lm", "logistic", "brlogistic", "poisson", 
                     "coxph", "compet", "casecohort", "nestedcc", "weibull", 
                     "dtsam", "flexsurv")) == FALSE) {
    stop(paste("Substantive model type ", smtype, " not recognised.", 
               sep = ""))
  }
  if (smtype == "flexsurv") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    if (!(all(sort(unique(d)) == c(0, 1))) & !(all(unique(d) == 
                                                   1))) {
      stop("Event indicator for flexsurv must be coded 0/1 for censoring/event.")
    }
    if ((sum(is.na(d)) + sum(is.na(originaldata[, timeCol]))) > 
        0) {
      stop("Event indicator and time variables should not have NAs.")
    }
  }
  else if (smtype == "coxph") { ## modified
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(in.subco)]
    d <- originaldata[, dCol]
    s <- originaldata[, subcoCol]
    subco.weight <- originaldata$weights
    list.times <- sort(unique(originaldata[, timeCol][originaldata[, dCol] == 1]))
    
    if (!(all(sort(unique(d)) == c(0, 1))) & !(all(unique(d) == 
                                                   1))) {
      stop("Event indicator for coxph must be coded 0/1 for censoring/event.")
    }
    if ((sum(is.na(d)) + sum(is.na(originaldata[, timeCol]))) > 
        0) {
      stop("Event indicator and time variables should not have NAs.")
    }
    nullMod <- survival::coxph(survival::Surv(originaldata[, 
                                                           timeCol], originaldata[, dCol]) ~ 1, control = survival::coxph.control(timefix = FALSE))
    basehaz <- survival::basehaz(nullMod)
    H0indices <- match(originaldata[, timeCol], basehaz[, 
                                                        2])
    rm(nullMod)
  }
  else if (smtype == "weibull") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
  }
  else if (smtype == "dtsam") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    cutPoints <- 1:max(originaldata[, timeCol])
    nTimePoints <- length(cutPoints)
    if (!all(unique(originaldata[, timeCol]) == floor(unique(originaldata[, 
                                                                          timeCol])))) {
      stop("Your time variable must only take positive integer values.")
    }
    if (any(unique(originaldata[, timeCol]) <= 0)) {
      stop("Your time variable must only take positive integer values.")
    }
    if (extraArgs$timeEffects == "factor") {
      if (!identical(sort(unique(originaldata[, timeCol][originaldata[, 
                                                                      dCol] == 1])), as.numeric(cutPoints))) {
        stop("You cannot fit a dtsam model with factor time effects since there are some periods with no events. See documentation for timeEffects argument for parametric alternatives")
      }
    }
  }
  else if (smtype == "compet") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula[[1]])[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula[[1]])[[2]][[3]][[2]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    numCauses <- length(smformula)
    H0 <- vector("list", numCauses)
    H0indices <- vector("list", numCauses)
    outcomeModBeta <- vector("list", numCauses)
    linpred <- vector("list", numCauses)
    for (cause in 1:numCauses) {
      nullMod <- survival::coxph(as.formula(paste(strsplit(smformula[[cause]], 
                                                           "~")[[1]][1], "~1")), originaldata, control = survival::coxph.control(timefix = FALSE))
      basehaz <- survival::basehaz(nullMod)
      H0[[cause]] <- basehaz[, 1]
      H0indices[[cause]] <- match(originaldata[, timeCol], 
                                  basehaz[, 2])
      linpred[[cause]] <- as.formula(smformula[[cause]])
    }
    rm(nullMod)
  }
  else if (smtype == "casecohort") {
    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                           extraArgs$in.subco]
    subcoMembers <- which(originaldata[, subcoCol] == 1)
    subco.weight <- ifelse(originaldata[, subcoCol] == 1, 
                           1/extraArgs$sampfrac, 0)
    entertimeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                               toString(as.formula(smformula)[[2]][[2]])]
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[3]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[4]])]
    outcomeCol <- c(entertimeCol, timeCol, dCol)
    d <- originaldata[, dCol]
    list.times <- sort(unique(originaldata[, timeCol][originaldata[, 
                                                                   dCol] == 1]))
    smcfcsid <- 1:n
    smformula2 <- paste(smformula, "+cluster(smcfcsid)", 
                        sep = "")
  }
  else if (smtype == "nestedcc") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                          toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                       toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    setCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                         toString(extraArgs$set)]
    nriskCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                           toString(extraArgs$nrisk)]
    eventCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% 
                                           toString(extraArgs$event)]
    exp1 <- as.formula(paste(smformula))[[2]]
    exp2 <- as.formula(smformula)[[3]][[2]]
    smformula2 <- paste(deparse(exp1), "~", deparse(exp2, 
                                                    width.cutoff = 500L))
    d <- originaldata[, eventCol]
    noncases <- which(originaldata[, eventCol] == 0)
    num.sampriskset <- ave(rep(1, dim(originaldata)[1]), 
                           originaldata[, setCol], FUN = function(x) sum(x))
  }
  else {
    outcomeCol <- which(colnames(originaldata) == as.formula(smformula)[[2]])
  }
  if ((smtype == "logistic") | (smtype == "brlogistic")) {
    if (is.numeric(originaldata[, outcomeCol]) == FALSE) {
      stop("For logistic substantive models the outcome variable must be numeric 0/1.")
    }
    else {
      if (all.equal(sort(unique(originaldata[, outcomeCol])), 
                    c(0, 1)) == FALSE) {
        stop("For logistic substantive models the outcome variable must be coded 0/1.")
      }
    }
  }
  if (smtype == "compet") {
    smcovnames <- attr(terms(as.formula(smformula[[1]])), 
                       "term.labels")
    for (cause in 2:numCauses) {
      smcovnames <- c(smcovnames, attr(terms(as.formula(smformula[[cause]])), 
                                       "term.labels"))
    }
    smcovnames <- unique(smcovnames)
  }
  else if (smtype == "nestedcc") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")[-length(attr(terms(as.formula(smformula)), 
                                                                                 "term.labels"))]
  }
  else {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% 
                                        smcovnames]
  partialVars <- which((method == "norm") | (method == "latnorm") | 
                         (method == "logreg") | (method == "poisson") | (method == 
                                                                           "podds") | (method == "mlogit") | (method == "brlogreg"))
  if (length(partialVars) == 0) {
    if (((smtype == "flexsurv") & (extraArgs$imputeTimes == 
                                   TRUE)) == FALSE) 
      stop("You have not specified any valid imputation methods in the method argument.")
  }
  for (colnum in 1:ncol(originaldata)) {
    if (method[colnum] != "") {
      if (colnum %in% outcomeCol) {
        stop(paste("An imputation method has been specified for ", 
                   colnames(originaldata)[colnum], ". Elements of the method argument corresponding to the outcome variable(s) should be empty.", 
                   sep = ""))
      }
      else {
        if (sum(r[, colnum]) == n) {
          stop(paste("An imputation method has been specified for ", 
                     colnames(originaldata)[colnum], ", but it appears to be fully observed.", 
                     sep = ""))
        }
      }
    }
    else {
      if (sum(r[, colnum]) < n) {
        if (((colnum %in% outcomeCol) == FALSE) & (sum(errorProneMatrix[, 
                                                                        colnum]) == 0)) {
          stop(paste("Variable ", colnames(originaldata)[colnum], 
                     " does not have an imputation method specified, yet appears to have missing values.", 
                     sep = ""))
        }
      }
    }
  }
  if ("latnorm" %in% method) {
    if (is.null(errorProneMatrix) == TRUE) {
      stop("If you specify method latnorm you must specify the errorProneMatrix argument.")
    }
    else {
      if (identical(dim(errorProneMatrix), c(length(originaldata), 
                                             length(originaldata))) == FALSE) {
        stop("The errorProneMatrix should be a square matrix with number of rows equal to the number of variables in the dataset.")
      }
      if (identical(sort(unique(as.vector(errorProneMatrix))), 
                    c(0, 1)) == FALSE) {
        stop("The errorProneMatrix should only consist of 0s and 1s.")
      }
      for (varNum in 1:length(method)) {
        if (method[varNum] == "latnorm") {
          if (sum(errorProneMatrix[varNum, ]) < 2) {
            stop("Each latnorm variable must have two or more error prone measurements specified in the errorProneMatrix argument.")
          }
        }
      }
      if (sum(colSums(errorProneMatrix) > 1) > 0) {
        stop("Each error-prone measurement should be allocated to exactly one latnorm variable.")
      }
    }
  }
  if (is.null(errorProneMatrix) == FALSE) {
    if (("latnorm" %in% method) == FALSE) {
      stop("If you specify errorProneMatrix then at least one variable must be imputed using latnorm.")
    }
  }
  fullObsVars <- which((colSums(r) == n) & (colnames(originaldata) %in% 
                                              smcovnames))
  passiveVars <- which((method != "") & (method != "norm") & 
                         (method != "logreg") & (method != "poisson") & (method != 
                                                                           "podds") & (method != "mlogit") & (method != "latnorm") & 
                         (method != "brlogreg"))
  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol], 
                                            collapse = ",")))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars], 
                                          collapse = ",")))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars], 
                                                 collapse = ",")))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars], 
                                                               collapse = ",")))
  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }
  rjFailCount <- 0
  flexsurvFailCount <- 0
  for (imp in 1:m) {
    print(paste("Imputation ", imp))
    if (length(partialVars) > 0) {
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (method[targetCol] == "latnorm") {
          errorProneCols <- which(errorProneMatrix[targetCol, 
          ] == 1)
          for (measure in 1:length(errorProneCols)) {
            if (sum(r[, errorProneCols[measure]]) < 
                n) {
              imputations[[imp]][r[, errorProneCols[measure]] == 
                                   0, errorProneCols[measure]] <- sample(imputations[[imp]][r[, 
                                                                                              errorProneCols[measure]] == 1, errorProneCols[measure]], 
                                                                         size = sum(r[, errorProneCols[measure]] == 
                                                                                      0), replace = TRUE)
            }
          }
          imputations[[imp]][, targetCol] <- apply(imputations[[imp]][, 
                                                                      errorProneCols], 1, mean)
        }
        else {
          imputations[[imp]][r[, targetCol] == 0, targetCol] <- sample(imputations[[imp]][r[, 
                                                                                            targetCol] == 1, targetCol], size = sum(r[, 
                                                                                                                                      targetCol] == 0), replace = TRUE)
        }
      }
    }
    if ((smtype == "lm") | (smtype == "logistic") | (smtype == 
                                                     "brlogistic") | (smtype == "poisson")) {
      if (sum(r[, outcomeCol]) < n) {
        if (imp == 1) {
          print("Imputing missing outcomes using specified substantive model.")
        }
        imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                                method, passiveVars)
        imputationNeeded <- (1:n)[r[, outcomeCol] == 
                                    0]
        if (smtype == "lm") {
          ymod <- stats::lm(as.formula(smformula), imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), 
                                   imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded, outcomeCol] <- rnorm(length(imputationNeeded), 
                                                                    outmodxb[imputationNeeded], sigmasq^0.5)
        }
        else if ((smtype == "logistic") | (smtype == 
                                           "brlogistic")) {
          if (smtype == "logistic") {
            ymod <- glm(as.formula(smformula), family = "binomial", 
                        imputations[[imp]])
          }
          else {
            ymod <- glm(as.formula(smformula), family = "binomial", 
                        imputations[[imp]], method = brglm2::brglmFit)
          }
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), 
                                   imputations[[imp]]) %*% beta
          prob <- expit(outmodxb[imputationNeeded])
          imputations[[imp]][imputationNeeded, outcomeCol] <- rbinom(length(imputationNeeded), 
                                                                     1, prob)
        }
        else if (smtype == "poisson") {
          ymod <- glm(as.formula(smformula), family = "poisson", 
                      imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), 
                                   imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded, outcomeCol] <- rpois(length(imputationNeeded), 
                                                                    exp(outmodxb[imputationNeeded]))
        }
      }
    }
    for (cyclenum in 1:numit) {
      if (noisy == TRUE) {
        print(paste("Iteration ", cyclenum))
      }
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                              method, passiveVars)
      if (length(partialVars) > 0) {
        for (var in 1:length(partialVars)) {
          targetCol <- partialVars[var]
          if (is.null(predictorMatrix)) {
            predictorCols <- c(partialVars[!partialVars %in% 
                                             targetCol], fullObsVars)
          }
          else {
            predictorCols <- which(predictorMatrix[targetCol, 
            ] == 1)
            predictorCols <- predictorCols[!predictorCols %in% 
                                             outcomeCol]
          }
          if ((imp == 1) & (cyclenum == 1)) {
            if (method[targetCol] == "latnorm") {
              print(paste("Imputing: ", colnames(imputations[[imp]])[targetCol], 
                          " using ", paste(colnames(imputations[[imp]])[c(predictorCols, 
                                                                          which(errorProneMatrix[targetCol, 
                                                                          ] == 1))], collapse = ","), " plus outcome", 
                          collapse = ","))
            }
            else {
              print(paste("Imputing: ", colnames(imputations[[imp]])[targetCol], 
                          " using ", paste(colnames(imputations[[imp]])[predictorCols], 
                                           collapse = ","), " plus outcome", 
                          collapse = ","))
            }
          }
          if (length(predictorCols) > 0) {
            xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], 
                                            "~", paste(colnames(imputations[[imp]])[predictorCols], 
                                                       collapse = "+"), sep = ""))
          }
          else {
            xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], 
                                            "~1", sep = ""))
          }
          if (smtype == "casecohort") {
            xmoddata <- imputations[[imp]][subcoMembers, 
            ]
          }
          else if (smtype == "nestedcc") {
            xmoddata <- imputations[[imp]][noncases, 
            ]
          }
          ## modified
          else if (smtype == "coxph") {
            xmoddata <- imputations[[imp]][which(s == 1), ]
          }
          else {
            xmoddata <- imputations[[imp]]
          }
          if (method[targetCol] == "norm") {
            xmod <- lm(xmodformula, data = xmoddata)
            beta <- xmod$coef
            sigmasq <- summary(xmod)$sigma^2
            newsigmasq <- (sigmasq * xmod$df)/rchisq(1, 
                                                     xmod$df)
            covariance <- (newsigmasq/sigmasq) * vcov(xmod)
            newbeta <- beta + MASS::mvrnorm(1, mu = rep(0, 
                                                        ncol(covariance)), Sigma = covariance)
            if ((smtype == "casecohort") | (smtype == "nestedcc") | 
                (smtype == "coxph")) { ## modified
              xfitted <- model.matrix(xmodformula, data = imputations[[imp]]) %*% 
                newbeta
            }
            else {
              xfitted <- model.matrix(xmod) %*% newbeta
            }
          }
          else if (method[targetCol] == "latnorm") {
            xmod <- lm(xmodformula, data = xmoddata)
            beta <- xmod$coef
            sigmasq <- summary(xmod)$sigma^2
            newsigmasq <- 1/rgamma(1, shape = ((n + 
                                                  1)/2), rate = ((n * sigmasq + 1)/2))
            covariance <- (newsigmasq/sigmasq) * vcov(xmod)
            newbeta <- beta + MASS::mvrnorm(1, mu = rep(0, 
                                                        ncol(covariance)), Sigma = covariance)
            if ((smtype == "casecohort") | (smtype == 
                                            "nestedcc")) {
              xfitted <- model.matrix(xmodformula, data = imputations[[imp]]) %*% 
                newbeta
            }
            else {
              xfitted <- model.matrix(xmod) %*% newbeta
            }
            errorProneCols <- which(errorProneMatrix[targetCol, 
            ] == 1)
            xmat <- matrix(imputations[[imp]][, targetCol], 
                           nrow = nrow(imputations[[imp]]), ncol = length(errorProneCols))
            uVec <- c(as.matrix(imputations[[imp]][, 
                                                   errorProneCols] - xmat))
            sigmausq <- mean(uVec^2)
            sum_ni <- n * length(errorProneCols)
            sigmausq <- 1/rgamma(1, shape = ((sum_ni + 
                                                1)/2), rate = ((sum_ni * sigmausq + 1)/2))
            for (measure in 1:length(errorProneCols)) {
              nToImpute <- n - sum(r[, errorProneCols[measure]])
              if (nToImpute > 0) {
                imputations[[imp]][r[, errorProneCols[measure]] == 
                                     0, errorProneCols[measure]] <- imputations[[imp]][r[, 
                                                                                         errorProneCols[measure]] == 0, targetCol] + 
                  rnorm(nToImpute, 0, sd = sqrt(sigmausq))
              }
            }
            wmean <- rowMeans(imputations[[imp]][, errorProneCols])
            lambda <- newsigmasq/(newsigmasq + sigmausq/length(errorProneCols))
            xfitted <- xfitted + lambda * (wmean - xfitted)
            newsigmasq <- rep(newsigmasq * (1 - lambda), 
                              n)
          }
          else if (method[targetCol] == "logreg") {
            xmod <- glm(xmodformula, family = "binomial", 
                        data = xmoddata)
            newbeta <- modPostDraw(xmod)
            if ((smtype == "casecohort") | (smtype == "nestedcc") | (smtype == "coxph")) {
              xfitted <- expit(model.matrix(xmodformula, 
                                            data = imputations[[imp]]) %*% newbeta)
            }
            else {
              xfitted <- expit(model.matrix(xmod) %*% 
                                 newbeta)
            }
          }
          else if (method[targetCol] == "brlogreg") {
            xmod <- glm(xmodformula, family = "binomial", 
                        data = xmoddata, method = brglm2::brglmFit)
            newbeta <- modPostDraw(xmod)
            if ((smtype == "casecohort") | (smtype == "nestedcc") | 
                (smtype == "coxph")) {
              xfitted <- expit(model.matrix(xmodformula, 
                                            data = imputations[[imp]]) %*% newbeta)
            }
            else {
              xfitted <- expit(model.matrix(xmod) %*% 
                                 newbeta)
            }
          }
          else if (method[targetCol] == "poisson") {
            xmod <- glm(xmodformula, family = "poisson", 
                        data = xmoddata)
            newbeta <- modPostDraw(xmod)
            if ((smtype == "casecohort") | (smtype == 
                                            "nestedcc")) {
              xfitted <- exp(model.matrix(xmodformula, 
                                          data = imputations[[imp]]) %*% newbeta)
            }
            else {
              xfitted <- exp(model.matrix(xmod) %*% 
                               newbeta)
            }
          }
          else if (method[targetCol] == "podds") {
            if (is.ordered(imputations[[imp]][, targetCol]) == 
                FALSE) 
              stop("Variables to be imputed using method podds must be stored as ordered factors.")
            xmod <- VGAM::vglm(xmodformula, VGAM::propodds, 
                               data = xmoddata)
            xmod.dummy <- VGAM::vglm(xmodformula, VGAM::propodds, 
                                     data = imputations[[imp]])
            newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, 
                                                        mu = rep(0, ncol(VGAM::vcov(xmod))), Sigma = VGAM::vcov(xmod))
            linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% 
                                 newbeta, byrow = TRUE, ncol = (nlevels(imputations[[imp]][, 
                                                                                           targetCol]) - 1))
            cumprobs <- cbind(1/(1 + exp(linpreds)), 
                              rep(1, nrow(linpreds)))
            xfitted <- cbind(cumprobs[, 1], cumprobs[, 
                                                     2:ncol(cumprobs)] - cumprobs[, 1:(ncol(cumprobs) - 
                                                                                         1)])
          }
          else if (method[targetCol] == "mlogit") {
            if (is.factor(imputations[[imp]][, targetCol]) == 
                FALSE) 
              stop("Variables to be imputed using method modds must be stored as factors.")
            xmod <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel = 1), 
                               data = xmoddata)
            xmod.dummy <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel = 1), 
                                     data = imputations[[imp]])
            newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, 
                                                        mu = rep(0, ncol(VGAM::vcov(xmod))), Sigma = VGAM::vcov(xmod))
            linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% 
                                 newbeta, byrow = TRUE, ncol = (nlevels(imputations[[imp]][, 
                                                                                           targetCol]) - 1))
            denom <- 1 + rowSums(exp(linpreds))
            xfitted <- cbind(1/denom, exp(linpreds)/denom)
          }
          if (noisy == TRUE) {
            print(summary(xmod))
          }
          if (smtype == "lm") {
            ymod <- lm(as.formula(smformula), imputations[[imp]])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq * ymod$df)/rchisq(1, 
                                                           ymod$df)
            covariance <- (outcomeModResVar/sigmasq) * 
              vcov(ymod)
            outcomeModBeta <- beta + MASS::mvrnorm(1, 
                                                   mu = rep(0, ncol(covariance)), Sigma = covariance)
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          else if ((smtype == "logistic") | (smtype == 
                                             "brlogistic")) {
            if (smtype == "logistic") {
              ymod <- glm(as.formula(smformula), family = "binomial", 
                          imputations[[imp]])
            }
            else {
              ymod <- glm(as.formula(smformula), family = "binomial", 
                          imputations[[imp]], method = brglm2::brglmFit)
            }
            outcomeModBeta <- modPostDraw(ymod)
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          else if (smtype == "dtsam") {
            longData <- survival::survSplit(as.formula(smformula), 
                                            data = imputations[[imp]], cut = cutPoints)
            if (extraArgs$timeEffects == "factor") {
              dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], 
                                    "~-1+factor(tstart)+", strsplit(smformula, 
                                                                    "~")[[1]][2], sep = "")
            }
            else if (extraArgs$timeEffects == "linear") {
              dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], 
                                    "~tstart+", strsplit(smformula, "~")[[1]][2], 
                                    sep = "")
            }
            else {
              dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], 
                                    "~tstart+I(tstart^2)+", strsplit(smformula, 
                                                                     "~")[[1]][2], sep = "")
            }
            ymod <- glm(as.formula(dtsamFormula), family = "binomial", 
                        data = longData)
            outcomeModBeta <- modPostDraw(ymod)
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          else if (smtype == "poisson") {
            ymod <- glm(as.formula(smformula), family = "poisson", 
                        imputations[[imp]])
            outcomeModBeta <- modPostDraw(ymod)
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          } ## modified
          else if (smtype == "coxph") {
            if (super == F) {
              ymod <- survival::coxph(as.formula(smformula), 
                                      data = imputations[[imp]], 
                                      control = survival::coxph.control(timefix = FALSE),
                                      cluster = imputations[[imp]]$id)
              outcomeModBeta <- modPostDraw(ymod)
              ymod$coefficients <- outcomeModBeta
              basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
              H0 <- basehaz[H0indices]
            }
            else if (super == T){
              ymod <- survival::coxph(as.formula(smformula), 
                                      data = imputations[[imp]], 
                                      control = survival::coxph.control(timefix = FALSE),
                                      cluster = imputations[[imp]]$id, weights = imputations[[imp]]$weights)
              outcomeModBeta <- modPostDraw(ymod)
              cumhaz.denom.elements <- exp(model.matrix(as.formula(smformula),imputations[[imp]])[,-1] %*% outcomeModBeta)
              cumhaz.denom <- sapply(list.times,function(x){
                sum(cumhaz.denom.elements[which(originaldata[,timeCol]>=x)] * subco.weight[which(originaldata[,timeCol]>=x)])
              })
              exp.func.denom <- cumsum(1/cumhaz.denom)
              H0.fun <- stepfun(list.times,c(0,exp.func.denom))
              H0 <- H0.fun(originaldata[,timeCol])
            }
            if (noisy==TRUE) {
              print(summary(ymod))
            }
          }
          else if (smtype == "flexsurv") {
            if (cyclenum == 1) {
              tryCatch({
                ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                                 imputations[[imp]], k = extraArgs$k, 
                                                 scale = "hazard")
              }, error = function(e) {
                flexsurvFailCount <<- flexsurvFailCount + 
                  1
              })
            }
            else {
              tryCatch({
                if (extraArgs$originalKnots == FALSE) {
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                                   imputations[[imp]], k = extraArgs$k, 
                                                   scale = "hazard", inits = flexsurvEsts)
                }
                else {
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                                   imputations[[imp]], scale = "hazard", 
                                                   knots = ymod$knots[2:(length(ymod$knots) - 
                                                                           1)], inits = flexsurvEsts)
                }
              }, error = function(e) {
                flexsurvFailCount <<- flexsurvFailCount + 
                  1
              })
            }
            flexsurvEsts <- ymod$res.t[, 1]
            outcomeModBeta <- as.numeric(flexsurv::normboot.flexsurvreg(ymod, 
                                                                        B = 1, raw = TRUE))
            ymod$res.t[, 1] <- outcomeModBeta
            if (noisy == TRUE) {
              print(ymod)
            }
          }
          else if (smtype == "weibull") {
            ymod <- survival::survreg(as.formula(smformula), 
                                      data = imputations[[imp]], dist = "weibull")
            outcomeModBeta <- c(coef(ymod), log(ymod$scale)) + 
              MASS::mvrnorm(1, mu = rep(0, ncol(vcov(ymod))), 
                            Sigma = vcov(ymod))
            weibullScale <- exp(utils::tail(outcomeModBeta, 
                                            1))
            outcomeModBeta <- utils::head(outcomeModBeta, 
                                          -1)
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          else if (smtype == "compet") {
            for (cause in 1:numCauses) {
              ymod <- survival::coxph(as.formula(smformula[[cause]]), 
                                      imputations[[imp]], control = survival::coxph.control(timefix = FALSE))
              outcomeModBeta[[cause]] <- modPostDraw(ymod)
              ymod$coefficients <- outcomeModBeta[[cause]]
              basehaz <- survival::basehaz(ymod, centered = FALSE)[, 
                                                                   1]
              H0[[cause]] <- basehaz[H0indices[[cause]]]
              if (noisy == TRUE) {
                print(summary(ymod))
              }
            }
          }
          else if (smtype == "casecohort") {
            ymod <- survival::coxph(as.formula(smformula2), 
                                    imputations[[imp]], control = survival::coxph.control(timefix = FALSE))
            outcomeModBeta <- modPostDraw(ymod)
            cumhaz.denom.elements <- exp(model.matrix(as.formula(smformula), 
                                                      imputations[[imp]])[, -1] %*% outcomeModBeta)
            cumhaz.denom <- sapply(list.times, function(x) {
              sum(cumhaz.denom.elements[which(originaldata[, 
                                                           timeCol] >= x)] * subco.weight[which(originaldata[, 
                                                                                                             timeCol] >= x)])
            })
            exp.func.denom <- cumsum(1/cumhaz.denom)
            H0.fun <- stepfun(list.times, c(0, exp.func.denom))
            H0 <- H0.fun(originaldata[, timeCol])
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          else if (smtype == "nestedcc") {
            ymod <- survival::coxph(as.formula(smformula), 
                                    imputations[[imp]], control = survival::coxph.control(timefix = FALSE))
            outcomeModBeta <- modPostDraw(ymod)
            explan.matrix <- model.matrix(ymod)
            cumbasehaz.denom <- exp(matrix(outcomeModBeta, 
                                           nrow = 1) %*% t(explan.matrix)) * originaldata[, 
                                                                                          nriskCol]/num.sampriskset
            cumbasehaz.denom <- ave(cumbasehaz.denom, 
                                    originaldata[, setCol], FUN = sum)[originaldata[, 
                                                                                    dCol] == 1]
            cumbasehaz.t <- originaldata[, timeCol][originaldata[, 
                                                                 dCol] == 1]
            H0 <- unlist(lapply(originaldata[, timeCol], 
                                function(x) {
                                  sum((1/cumbasehaz.denom)[cumbasehaz.t <= 
                                                             x])
                                }))
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
          if ((imp == 1) & (cyclenum == 1) & (var == 
                                              1)) {
            if (smtype == "compet") {
              totalCoefVec <- outcomeModBeta[[1]]
              for (cause in 2:numCauses) {
                totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
              }
              smCoefIter <- array(0, dim = c(m, length(totalCoefVec), 
                                             numit))
            }
            else {
              smCoefIter <- array(0, dim = c(m, length(outcomeModBeta), 
                                             numit))
            }
          }
          if (var == length(partialVars)) {
            if (smtype == "compet") {
              totalCoefVec <- outcomeModBeta[[1]]
              for (cause in 2:numCauses) {
                totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
              }
              smCoefIter[imp, , cyclenum] <- totalCoefVec
            }
            else {
              smCoefIter[imp, , cyclenum] <- outcomeModBeta
            }
          }
          imputationNeeded <- (1:n)[r[, targetCol] == 
                                      0]
          if ((method[targetCol] == "logreg") | (method[targetCol] == 
                                                 "podds") | (method[targetCol] == "mlogit") | 
              (method[targetCol] == "brlogreg")) {
            if ((method[targetCol] == "logreg") | (method[targetCol] == 
                                                   "brlogreg")) {
              numberOutcomes <- 2
              fittedMean <- cbind(1 - xfitted, xfitted)
            }
            else {
              numberOutcomes <- nlevels(imputations[[imp]][, 
                                                           targetCol])
              fittedMean <- xfitted
            }
            outcomeDensCovDens <- array(dim = c(length(imputationNeeded), 
                                                numberOutcomes), 0)
            for (xMisVal in 1:numberOutcomes) {
              if ((method[targetCol] == "logreg") | 
                  (method[targetCol] == "brlogreg")) {
                if (is.factor(imputations[[imp]][, targetCol]) == 
                    TRUE) {
                  valToImpute <- levels(imputations[[imp]][, 
                                                           targetCol])[xMisVal]
                }
                else {
                  valToImpute <- xMisVal - 1
                }
              }
              else {
                valToImpute <- levels(imputations[[imp]][, 
                                                         targetCol])[xMisVal]
              }
              imputations[[imp]][imputationNeeded, targetCol] <- valToImpute
              imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                                      method, passiveVars)
              if (smtype == "lm") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                deviation <- imputations[[imp]][imputationNeeded, 
                                                outcomeCol] - outmodxb[imputationNeeded]
                outcomeDens <- dnorm(deviation, mean = 0, 
                                     sd = outcomeModResVar^0.5)
              }
              else if ((smtype == "logistic") | (smtype == 
                                                 "brlogistic")) {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                prob <- expit(outmodxb[imputationNeeded])
                outcomeDens <- prob * imputations[[imp]][imputationNeeded, 
                                                         outcomeCol] + (1 - prob) * (1 - imputations[[imp]][imputationNeeded, 
                                                                                                            outcomeCol])
              }
              else if (smtype == "dtsam") {
                outcomeDens <- dtsamOutcomeDens(imputations[[imp]], 
                                                extraArgs$timeEffects, outcomeModBeta, 
                                                nTimePoints, smformula, timeCol, dCol)[imputationNeeded]
              }
              else if (smtype == "poisson") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                outcomeDens <- dpois(imputations[[imp]][imputationNeeded, 
                                                        outcomeCol], exp(outmodxb[imputationNeeded]))
              }
              else if (smtype == "weibull") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                outcomeDens <- (1 - d[imputationNeeded]) * 
                  (1 - survival::psurvreg(imputations[[imp]][imputationNeeded, 
                                                             timeCol], mean = outmodxb[imputationNeeded], 
                                          scale = weibullScale)) + d[imputationNeeded] * 
                  (survival::dsurvreg(imputations[[imp]][imputationNeeded, 
                                                         timeCol], mean = outmodxb[imputationNeeded], 
                                      scale = weibullScale))
              }
              else if ((smtype == "coxph") | (smtype == 
                                              "casecohort")) { ## modified
                outmodxb <- model.matrix(as.formula(smformula),
                                         imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                outcomeDens <- exp(-H0[imputationNeeded] * 
                                     exp(outmodxb[imputationNeeded])) * 
                  (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
              }
              else if (smtype == "flexsurv") {
                survEst <- summary(ymod, newdata = imputations[[imp]][imputationNeeded, 
                ], type = "survival", ci = FALSE, 
                t = imputations[[imp]][imputationNeeded, 
                                       timeCol], cross = FALSE, tidy = TRUE)
                survEst <- as.matrix(survEst)[, "est"]
                hazEst <- summary(ymod, newdata = imputations[[imp]][imputationNeeded, 
                ], type = "hazard", ci = FALSE, t = imputations[[imp]][imputationNeeded, 
                                                                       timeCol], cross = FALSE, tidy = TRUE)
                hazEst <- as.matrix(hazEst)[, "est"]
                outcomeDens <- survEst * (hazEst^imputations[[imp]][imputationNeeded, 
                                                                    dCol])
              }
              else if (smtype == "nestedcc") {
                outmodxb <- model.matrix(as.formula(smformula2), 
                                         imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                outcomeDens <- exp(-H0[imputationNeeded] * 
                                     exp(outmodxb[imputationNeeded])) * 
                  (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
              }
              else if (smtype == "compet") {
                outcomeDens <- rep(1, length(imputationNeeded))
                for (cause in 1:numCauses) {
                  outmodxb <- model.matrix(linpred[[cause]], 
                                           imputations[[imp]])
                  outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                    as.matrix(outcomeModBeta[[cause]])
                  outcomeDens <- outcomeDens * exp(-H0[[cause]][imputationNeeded] * 
                                                     exp(outmodxb[imputationNeeded])) * 
                    (exp(outmodxb[imputationNeeded])^(d[imputationNeeded] == 
                                                        cause))
                }
              }
              outcomeDensCovDens[, xMisVal] <- outcomeDens * 
                fittedMean[imputationNeeded, xMisVal]
            }
            directImpProbs <- outcomeDensCovDens/rowSums(outcomeDensCovDens)
            if ((method[targetCol] == "logreg") | (method[targetCol] == 
                                                   "brlogreg")) {
              directImpProbs <- directImpProbs[, 2]
              if (is.factor(imputations[[imp]][, targetCol]) == 
                  TRUE) {
                imputations[[imp]][imputationNeeded, 
                                   targetCol] <- levels(imputations[[imp]][, 
                                                                           targetCol])[1]
                imputations[[imp]][imputationNeeded, 
                                   targetCol][rbinom(length(imputationNeeded), 
                                                     1, directImpProbs) == 1] <- levels(imputations[[imp]][, 
                                                                                                           targetCol])[2]
              }
              else {
                imputations[[imp]][imputationNeeded, 
                                   targetCol] <- rbinom(length(imputationNeeded), 
                                                        1, directImpProbs)
              }
            }
            else {
              imputations[[imp]][imputationNeeded, targetCol] <- levels(imputations[[imp]][, 
                                                                                           targetCol])[apply(directImpProbs, 1, 
                                                                                                             catdraw)]
            }
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                                    method, passiveVars)
          }
          else {
            firstTryLimit <- 25
            j <- 1
            while ((length(imputationNeeded) > 0) & 
                   (j < firstTryLimit)) {
              if ((method[targetCol] == "norm") | (method[targetCol] == 
                                                   "latnorm")) {
                imputations[[imp]][imputationNeeded, 
                                   targetCol] <- rnorm(length(imputationNeeded), 
                                                       xfitted[imputationNeeded], newsigmasq^0.5)
              }
              else if (method[targetCol] == "poisson") {
                imputations[[imp]][imputationNeeded, 
                                   targetCol] <- rpois(length(imputationNeeded), 
                                                       xfitted[imputationNeeded])
              }
              imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                                      method, passiveVars)
              uDraw <- runif(length(imputationNeeded))
              if (smtype == "lm") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                deviation <- imputations[[imp]][imputationNeeded, 
                                                outcomeCol] - outmodxb[imputationNeeded]
                reject <- 1 * (log(uDraw) > -(deviation^2)/(2 * 
                                                              array(outcomeModResVar, dim = c(length(imputationNeeded), 
                                                                                              1))))
              }
              else if ((smtype == "logistic") | (smtype == 
                                                 "brlogistic")) {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                prob <- expit(outmodxb[imputationNeeded])
                prob <- prob * imputations[[imp]][imputationNeeded, 
                                                  outcomeCol] + (1 - prob) * (1 - imputations[[imp]][imputationNeeded, 
                                                                                                     outcomeCol])
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "dtsam") {
                prob <- dtsamOutcomeDens(imputations[[imp]], 
                                         extraArgs$timeEffects, outcomeModBeta, 
                                         nTimePoints, smformula, timeCol, dCol)[imputationNeeded]
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "poisson") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                prob <- dpois(imputations[[imp]][imputationNeeded, 
                                                 outcomeCol], exp(outmodxb[imputationNeeded]))
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "weibull") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]]) %*% outcomeModBeta
                s_t <- 1 - survival::psurvreg(imputations[[imp]][imputationNeeded, 
                                                                 timeCol], mean = outmodxb[imputationNeeded], 
                                              scale = weibullScale)
                prob <- -exp(1) * log(s_t) * s_t
                prob <- d[imputationNeeded] * prob + 
                  (1 - d[imputationNeeded]) * s_t
                reject <- 1 * (uDraw > prob)
              }
              else if ((smtype == "coxph") | (smtype == 
                                              "casecohort")) { ## modified
                outmodxb <- model.matrix(as.formula(smformula), 
                                         imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                s_t <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))
                prob <- exp(1 + outmodxb[imputationNeeded] - 
                              (H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))) * 
                  H0[imputationNeeded]
                prob <- d[imputationNeeded] * prob + 
                  (1 - d[imputationNeeded]) * s_t
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "flexsurv") {
                survEst <- summary(ymod, newdata = imputations[[imp]][imputationNeeded, 
                ], type = "survival", ci = FALSE, 
                t = imputations[[imp]][imputationNeeded, 
                                       timeCol], cross = FALSE, tidy = TRUE)
                survEst <- as.matrix(survEst)[, "est"]
                hazEst <- summary(ymod, newdata = imputations[[imp]][imputationNeeded, 
                ], type = "hazard", ci = FALSE, t = imputations[[imp]][imputationNeeded, 
                                                                       timeCol], cross = FALSE, tidy = TRUE)
                hazEst <- as.matrix(hazEst)[, "est"]
                cumhazEst <- summary(ymod, newdata = imputations[[imp]][imputationNeeded, 
                ], type = "cumhaz", ci = FALSE, t = imputations[[imp]][imputationNeeded, 
                                                                       timeCol], cross = FALSE, tidy = TRUE)
                cumhazEst <- as.matrix(cumhazEst)[, 
                                                  "est"]
                prob <- imputations[[imp]][imputationNeeded, 
                                           dCol] * (survEst * exp(1) * cumhazEst) + 
                  (1 - imputations[[imp]][imputationNeeded, 
                                          dCol]) * survEst
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "nestedcc") {
                outmodxb <- model.matrix(as.formula(smformula2), 
                                         imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                s_t <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))
                prob <- exp(1 + outmodxb[imputationNeeded] - 
                              (H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))) * 
                  H0[imputationNeeded]
                prob <- d[imputationNeeded] * prob + 
                  (1 - d[imputationNeeded]) * s_t
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "compet") {
                prob <- rep(1, length(imputationNeeded))
                for (cause in 1:numCauses) {
                  outmodxb <- model.matrix(linpred[[cause]], 
                                           imputations[[imp]])
                  outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                    as.matrix(outcomeModBeta[[cause]])
                  prob <- prob * exp(-H0[[cause]][imputationNeeded] * 
                                       exp(outmodxb[imputationNeeded])) * 
                    (H0[[cause]][imputationNeeded] * 
                       exp(1 + outmodxb[imputationNeeded]))^(d[imputationNeeded] == 
                                                               cause)
                }
                reject <- 1 * (uDraw > prob)
              }
              imputationNeeded <- imputationNeeded[reject == 
                                                     1]
              j <- j + 1
            }
            for (i in imputationNeeded) {
              tempData <- imputations[[imp]][i, ]
              tempData <- tempData[rep(1, rjlimit), 
              ]
              if (method[targetCol] == "norm") {
                tempData[, targetCol] <- rnorm(rjlimit, 
                                               xfitted[i], newsigmasq^0.5)
              }
              else if ((method[targetCol] == "logreg") | 
                       (method[targetCol] == "brlogreg")) {
                tempData[, targetCol] <- rbinom(rjlimit, 
                                                size = 1, xfitted[i])
              }
              else if (method[targetCol] == "poisson") {
                tempData[, targetCol] <- rpois(rjlimit, 
                                               xfitted[i])
              }
              else if (method[targetCol] == "latnorm") {
                tempData[, targetCol] <- rnorm(rjlimit, 
                                               xfitted[i], newsigmasq[i]^0.5)
              }
              tempData <- updatePassiveVars(tempData, 
                                            method, passiveVars)
              uDraw <- runif(rjlimit)
              if (smtype == "lm") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         tempData) %*% outcomeModBeta
                deviation <- tempData[, outcomeCol] - 
                  outmodxb
                reject <- 1 * (log(uDraw) > -(deviation^2)/(2 * 
                                                              array(outcomeModResVar, dim = c(rjlimit, 
                                                                                              1))))
              }
              else if ((smtype == "logistic") | (smtype == 
                                                 "brlogistic")) {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         tempData) %*% outcomeModBeta
                prob <- expit(outmodxb)
                prob <- prob * tempData[, outcomeCol] + 
                  (1 - prob) * (1 - tempData[, outcomeCol])
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "dtsam") {
                prob <- dtsamOutcomeDens(tempData, extraArgs$timeEffects, 
                                         outcomeModBeta, nTimePoints, smformula, 
                                         timeCol, dCol)
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "poisson") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         tempData) %*% outcomeModBeta
                prob <- dpois(tempData[, outcomeCol], 
                              exp(outmodxb))
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "weibull") {
                outmodxb <- model.matrix(as.formula(smformula), 
                                         tempData) %*% outcomeModBeta
                s_t <- 1 - survival::psurvreg(tempData[, 
                                                       timeCol], mean = outmodxb, scale = weibullScale)
                if (d[i] == 1) {
                  prob <- -exp(1) * log(s_t) * s_t
                  prob[is.na(prob)] <- 0
                }
                else {
                  prob <- s_t
                }
                reject <- 1 * (uDraw > prob)
              }
              else if ((smtype == "coxph") | (smtype == 
                                              "casecohort")) { ## modified
                outmodxb <- model.matrix(as.formula(smformula), 
                                         tempData)
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                s_t <- exp(-H0[i] * exp(outmodxb))
                prob <- exp(1 + outmodxb - (H0[i] * 
                                              exp(outmodxb))) * H0[i]
                prob <- d[i] * prob + (1 - d[i]) * s_t
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "flexsurv") {
                survEst <- summary(ymod, newdata = tempData, 
                                   type = "survival", ci = FALSE, t = tempData[, 
                                                                               timeCol], cross = FALSE, tidy = TRUE)
                survEst <- as.matrix(survEst)[, "est"]
                hazEst <- summary(ymod, newdata = tempData, 
                                  type = "hazard", ci = FALSE, t = tempData[, 
                                                                            timeCol], cross = FALSE, tidy = TRUE)
                hazEst <- as.matrix(hazEst)[, "est"]
                cumhazEst <- summary(ymod, newdata = tempData, 
                                     type = "cumhaz", ci = FALSE, t = tempData[, 
                                                                               timeCol], cross = FALSE, tidy = TRUE)
                cumhazEst <- as.matrix(cumhazEst)[, 
                                                  "est"]
                prob <- imputations[[imp]][i, dCol] * 
                  (survEst * exp(1) * cumhazEst) + (1 - 
                                                      imputations[[imp]][i, dCol]) * survEst
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "nestedcc") {
                outmodxb <- model.matrix(as.formula(smformula2), 
                                         tempData)
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                  as.matrix(outcomeModBeta)
                s_t <- exp(-H0[i] * exp(outmodxb))
                prob <- exp(1 + outmodxb - (H0[i] * 
                                              exp(outmodxb))) * H0[i]
                prob <- d[i] * prob + (1 - d[i]) * s_t
                reject <- 1 * (uDraw > prob)
              }
              else if (smtype == "compet") {
                prob <- rep(1, rjlimit)
                for (cause in 1:numCauses) {
                  outmodxb <- model.matrix(linpred[[cause]], 
                                           tempData)
                  outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% 
                    as.matrix(outcomeModBeta[[cause]])
                  prob <- prob * exp(-H0[[cause]][i] * 
                                       exp(outmodxb)) * (H0[[cause]][i] * 
                                                           exp(1 + outmodxb))^(d[i] == cause)
                }
                reject <- 1 * (uDraw > prob)
              }
              if (sum(reject) < rjlimit) {
                imputations[[imp]][i, targetCol] <- tempData[reject == 
                                                               0, targetCol][1]
              }
              else {
                rjFailCount <- rjFailCount + 1
              }
            }
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], 
                                                    method, passiveVars)
          }
        }
      }
      if ((smtype == "lm") | (smtype == "logistic") | 
          (smtype == "brlogistic")) {
        if (sum(r[, outcomeCol]) < n) {
          imputationNeeded <- (1:n)[r[, outcomeCol] == 
                                      0]
          if (smtype == "lm") {
            ymod <- lm(as.formula(smformula), imputations[[imp]][r[, 
                                                                   outcomeCol] == 1, ])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq * ymod$df)/rchisq(1, 
                                                           ymod$df)
            covariance <- (outcomeModResVar/sigmasq) * 
              vcov(ymod)
            outcomeModBeta <- beta + MASS::mvrnorm(1, 
                                                   mu = rep(0, ncol(covariance)), Sigma = covariance)
            outmodxb <- model.matrix(as.formula(smformula), 
                                     imputations[[imp]]) %*% outcomeModBeta
            imputations[[imp]][imputationNeeded, outcomeCol] <- rnorm(length(imputationNeeded), 
                                                                      outmodxb[imputationNeeded], sigmasq^0.5)
          }
          else if ((smtype == "logistic") | (smtype == 
                                             "brlogistic")) {
            if (smtype == "logistic") {
              ymod <- glm(as.formula(smformula), family = "binomial", 
                          imputations[[imp]][r[, outcomeCol] == 
                                               1, ])
            }
            else {
              ymod <- glm(as.formula(smformula), family = "binomial", 
                          imputations[[imp]][r[, outcomeCol] == 
                                               1, ], method = brglm2::brglmFit)
            }
            outcomeModBeta <- modPostDraw(ymod)
            outmodxb <- model.matrix(as.formula(smformula), 
                                     imputations[[imp]]) %*% outcomeModBeta
            prob <- expit(outmodxb[imputationNeeded])
            imputations[[imp]][imputationNeeded, outcomeCol] <- rbinom(length(imputationNeeded), 
                                                                       1, prob)
          }
        }
      }
      else if (smtype == "flexsurv") {
        if (extraArgs$imputeTimes == TRUE) {
          if (cyclenum == 1) {
            tryCatch({
              ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                               imputations[[imp]], k = extraArgs$k, 
                                               scale = "hazard")
            }, error = function(e) {
              flexsurvFailCount <<- flexsurvFailCount + 
                1
            })
          }
          else {
            tryCatch({
              if (extraArgs$originalKnots == FALSE) {
                ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                                 imputations[[imp]], k = extraArgs$k, 
                                                 scale = "hazard", inits = flexsurvEsts)
              }
              else {
                ymod <- flexsurv::flexsurvspline(as.formula(smformula), 
                                                 imputations[[imp]], scale = "hazard", 
                                                 knots = ymod$knots[2:(length(ymod$knots) - 
                                                                         1)], inits = flexsurvEsts)
              }
            }, error = function(e) {
              flexsurvFailCount <<- flexsurvFailCount + 
                1
            })
          }
          flexsurvEsts <- ymod$res.t[, 1]
          outcomeModBeta <- as.numeric(flexsurv::normboot.flexsurvreg(ymod, 
                                                                      B = 1, raw = TRUE))
          if (exists("smCoefIter") == FALSE) {
            smCoefIter <- array(0, dim = c(m, length(outcomeModBeta), 
                                           numit))
          }
          smCoefIter[imp, , cyclenum] <- outcomeModBeta
          ymod$res.t[, 1] <- outcomeModBeta
          if (noisy == TRUE) {
            print(ymod)
          }
          if (is.null(extraArgs$censtime)) {
            timeImp <- simulate(ymod, nsim = 1, newdata = imputations[[imp]][d == 
                                                                               0, ], start = originaldata[d == 0, timeCol])
            imputations[[imp]][d == 0, timeCol] <- timeImp$time_1
            imputations[[imp]][, dCol] <- 1
          }
          else {
            timeImp <- simulate(ymod, nsim = 1, newdata = imputations[[imp]][d == 
                                                                               0, ], start = originaldata[d == 0, timeCol], 
                                censtime = extraArgs$censtime)
            imputations[[imp]][d == 0, timeCol] <- timeImp$time_1
            imputations[[imp]][d == 0, dCol] <- timeImp$event_1
          }
        }
      }
    }
  }
  if (rjFailCount > 0) {
    warning(paste("Rejection sampling failed ", rjFailCount, 
                  " times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.", 
                  sep = ""))
  }
  if (flexsurvFailCount > 0) {
    warning(paste("Flexsurv fit failed ", flexsurvFailCount, 
                  " times (across all variables, iterations, and imputations). See documentation for more details.", 
                  sep = ""))
  }
  res <- list(impDatasets = imputations, smCoefIter = smCoefIter, 
              smInfo = list(smtype = smtype, smformula = smformula), 
              extraArgs = extraArgs)
  class(res) <- "smcfcs"
  return(res)
}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[, i] <- with(data, eval(parse(text = method[i])))
  }
  data
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu = rep(0, ncol(varcov)), Sigma = varcov)
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1, size = 1, prob = prob) == 1]
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}