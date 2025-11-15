######################################################################################
#### Functions for super sampling and weight calibration ####
####################################################################################
# Require full cohort data for `data` argument

supersampling <- function(data, sub.smformula, in.subco, n.super, random, design, stratum){
  if (missing(random)) {
    stop("Argument 'random' must be explicitly specified (TRUE or FALSE).")
  }
  if (!is.logical(random) || length(random) != 1) {
    stop("'random' must be a single logical value (TRUE or FALSE).")
  }
  ind.fail <- as.character(as.formula(sub.smformula)[[2]][[3]])
  idx_pool <- which(data[[in.subco]] == 0 & data[[ind.fail]] == 0) # Indices for super sampling
  idx_cc <- which(data[[in.subco]] == 1 | data[[ind.fail]] == 1) # Indices for the case-cohort sample
  idx_case <- which(data[[ind.fail]] == 1) # Indices for the cases
  idx_ctrl <- setdiff(idx_cc, idx_case) # Indices for the controls
  N <- nrow(data)
  # Weights
  wgt <- rep(0, N)
  wgt[idx_case] <- 1
  if (random == T){
    if (design == "case.cohort"){
      idx_ss <- sample(idx_pool, size = n.super) # Randomly sample the super sample
      data[[in.subco]][idx_ss] <- 1 # Update the indicator for subcohort by giving 1 to the super sample
      
      subco.cases <- which(data[[in.subco]]==1 & data[[ind.fail]]==1)
      idx_full <- c(idx_ss, idx_cc)
      idx_nc <- setdiff(idx_full, idx_case)
      den <- sum(data[[in.subco]])-length(subco.cases)
      stopifnot(den > 0L)
      wgt[idx_nc] <- (N-length(idx_case))/den
    }
    else if (design == "stratified.case.cohort") {
      if(!require(sampling)){
        install.packages("sampling")
        library(sampling)
      }
      pool_df <- data[idx_pool, , drop = FALSE]
      ord <- order(pool_df[[stratum]])
      pool_df <- pool_df[ord, , drop = FALSE]
      pool_index <- idx_pool[ord]
      
      lev <- levels(pool_df[[stratum]])
      if (!is.null(names(n.super))) {
        n.super <- as.integer(n.super[lev])
      } else if (length(n.super) != length(lev)) {
        stop("`n.super` length must equal the number of strata in the pool.")
      }
      
      valid_strata <- lev[n.super > 0]
      valid_sizes  <- n.super[n.super > 0]
      keep <- pool_df[[stratum]] %in% valid_strata
      pool_df_valid <- pool_df[keep, , drop = FALSE]
      pool_index_valid <- pool_index[keep]
      pool_df_valid[[stratum]] <- droplevels(pool_df_valid[[stratum]])
      
      lev_valid <- levels(pool_df_valid[[stratum]])
      size_vec  <- setNames(as.integer(valid_sizes), valid_strata)
      size_vec  <- size_vec[lev_valid]
      
      tab_valid <- table(pool_df_valid[[stratum]])
      if (any(size_vec > as.vector(tab_valid))) {
        stop("Requested `n.super` exceeds available units in at least one stratum.")
      }
      strata_sample <- sampling::strata(pool_df_valid, stratanames = stratum,
                                        size = size_vec, method = "srswor") # MUST specify "sampling::"strata
      idx_ss <- pool_index_valid[strata_sample$ID_unit]
      data[[in.subco]][idx_ss] <- 1 # Update the indicator for subcohort by giving 1 to the super sample
      
      subco.cases <- which(data[[in.subco]]==1 & data[[ind.fail]]==1)
      subco.cases <- split(subco.cases, data[[stratum]][subco.cases])
      cases <- split(idx_case, data[[stratum]][idx_case])
      subco <- which(data[[in.subco]]==1)
      subco <- split(subco, data[[stratum]][subco])
      N_str <- as.vector(table(data[[stratum]]))
      D_str <- lengths(cases)
      idx_full <- c(idx_ss, idx_cc)
      idx_nc <- setdiff(idx_full, idx_case)
      idx_nc_str <- split(idx_nc, data[[stratum]][idx_nc])
      wgt_vec <- (N_str-D_str)/(lengths(subco)-lengths(subco.cases))
      for (i in 1:length(idx_nc_str)) {
        wgt[idx_nc_str[[i]]] <- wgt_vec[i]
      }
    }
    wgt <- wgt[idx_full]
    return(structure(list(indices = idx_full, weights = wgt),
                     class = "supersampling"))
  }
  else {
    if (is.null(sub.smformula)) {
      stop("`sub.smformula` must be provided unless `random`==TRUE")
    }
    if(!require(sampling)){
      install.packages("sampling")
      library(sampling)
    }
    if(!require(BalancedSampling)){
      install.packages("BalancedSampling")
      library(BalancedSampling)
    }
    fit <- survival::coxph(as.formula(sub.smformula), control = survival::coxph.control(timefix = FALSE), data = data)
    dfbetas <- resid(fit, type="dfbetas")
    dfbetas <- as.matrix(dfbetas)
    abs.dfbetas <- abs(dfbetas)
    eps <- min(abs.dfbetas[abs.dfbetas>0]) * 1e-6
    db <- sqrt(rowSums(dfbetas^2)) + eps
    
    if (design == "case.cohort") {
      # Calculate inclusion probabilities
      incl_prob <- sampling::inclusionprobabilities(db[idx_pool], n.super)
      X <- cbind(incl_prob, dfbetas[idx_pool,])
      pos_ss <- BalancedSampling::cube(incl_prob, X)
      idx_ss <- idx_pool[pos_ss]
      # Indices of super sample and case-cohort sample
      idx_full <- c(idx_ss, idx_cc)
      
      # Weights
      nr <- nrow(data[data[[in.subco]]==1,])
      D <- length(idx_case)
      
      # Weights
      wgt[idx_ctrl] <- N/nr
      wgt[idx_ss] <- 1 / incl_prob[pos_ss]
      stopifnot(all(is.finite(wgt)))
      wgt <- wgt[idx_full]
      
      v_ctrl <- as.integer(seq(1,N) %in% idx_ctrl)
      v_ss <- as.integer(seq(1,N) %in% idx_ss)
      v_c <- as.integer(seq(1,N) %in% idx_case)
      
      db0 <- sum(db[idx_ctrl])
      db1 <- sum(db[idx_ss])
      Xs_cal <- cbind(cal1 = v_ctrl, cal2 = v_ss, cal3 = v_c)
      Xs_cal  <- Xs_cal[idx_full, , drop = FALSE]
      total_cal <- c((N-D)*db0/(db0+db1), (N-D)*db1/(db0+db1), D)
    }
    else if (design == "stratified.case.cohort") {
      n.pool <- table(data[[stratum]][idx_pool])
      idx_pool_str <- split(idx_pool, data[[stratum]][idx_pool])
      valid <- names(n.super)[n.super > 0 & n.pool > 0]
      ss_list <- lapply(valid, function(lev) {
        idxs <- idx_pool_str[[lev]]
        n_k  <- n.super[[lev]]
        n_p  <- n.pool[[lev]]
        
        if (n_k == 0) {
          return(list(idx = integer(0), wgt = numeric(0)))
        } else if (n_k >= n_p) {
          return(list(idx = idxs, wgt = setNames(rep(1, length(idxs)), idxs)))
        } else {
          pik <- sampling::inclusionprobabilities(db[idxs], n_k)
          X <- cbind(pik, dfbetas[idxs, , drop = FALSE])
          pos_ss <- BalancedSampling::cube(pik, X)
          return(list(idx = idxs[pos_ss], wgt = setNames(1 / pik[pos_ss], idxs[pos_ss]), p = pik))
        }
      })
      
      p_list <- lapply(ss_list, function(x) x$p)
      p <- unlist(p_list, use.names = FALSE)
      
      # Extract the indices for super sample
      idx_list <- lapply(ss_list, function(x) x$idx)
      idx_ss <- unlist(idx_list, use.names = FALSE)
      
      # Indices of super sample and case-cohort sample
      idx_full <- c(idx_ss, idx_cc)
      idx_nc <- setdiff(idx_full, idx_case)
      idx_nc_str <- split(idx_nc, data[[stratum]][idx_nc])
      N_str <- as.vector(table(data[[stratum]]))
      subco <- which(data[[in.subco]]==1)
      subco <- split(subco, data[[stratum]][subco])
      nr_str <- lengths(subco)
      D <- length(idx_case)
      
      wgt_vec2 <- N_str/nr_str
      for (i in 1:length(idx_nc_str)) {
        wgt[idx_nc_str[[i]]] <- wgt_vec2[i]
      }
      wgt_list <- lapply(ss_list, function(x) x$wgt)
      wgt_ss <- unlist(wgt_list, use.names = FALSE)
      wgt[idx_ss] <- wgt_ss
      wgt <- wgt[idx_full]
      
      v_ctrl <- as.integer(seq(1,N) %in% idx_ctrl)
      v_ss <- as.integer(seq(1,N) %in% idx_ss)
      v_c <- as.integer(seq(1,N) %in% idx_case)
      
      db0 <- sum(db[idx_ctrl])
      db1 <- sum(db[idx_ss])
      Xs_cal <- cbind(cal1 = v_ctrl, cal2 = v_ss, cal3 = v_c)
      Xs_cal  <- Xs_cal[idx_full, , drop = FALSE]
      total_cal <- c((N-D)*db0/(db0+db1), (N-D)*db1/(db0+db1), D)
    }
    g <- sampling::calib(
      Xs       = Xs_cal,
      d        = wgt,
      total    = total_cal,
      method   = "raking",
      max_iter = 500
    )
    return(structure(list(indices = idx_full, weights = wgt * g),
                     class = "super.sampling"))
  }
}