######################################################################################
#### Functions for super sampling and weight calibration ####
####################################################################################

add.ph2var <- function(fit, df, smformula, in.subco, cohort.size, stratum=NULL, cch.method, robust=NULL) {
  ind.fail <- as.character(as.formula(smformula)[[2]][[3]])
  subco.cases <- which(df[[in.subco]]==1 & df[[ind.fail]]==1)
  cases <- which(df[[ind.fail]]==1)
  subco.noncase <- which(df[[in.subco]]==1 & df[[ind.fail]]==0)
  subco <- which(df[[in.subco]]==1)
  if (cch.method == "LinYing"){
    wgt <- df$weights
    db <- as.matrix(resid(fit, type = "dfbeta"))
    db0 <- db[subco.noncase,,drop=F]
    dbm <- apply(db0,2,mean)
    db0 <- sweep(db0,2,dbm)
    
    fit$phase2var <- (1 - ((length(subco)-length(subco.cases))/(cohort.size-length(cases))))*crossprod(db0)
    fit$naive.var <- fit$naive.var + fit$phase2var
    
    if (robust)
      fit$var <- crossprod(db, db/wgt) + fit$phase2var
    else
      fit$var <- fit$naive.var
    return(fit)
  }
  else if (cch.method == "II.Borgan"){
    stratum <- as.numeric(df[[stratum]])
    nstrata <- max(stratum)
    score <- as.matrix(resid(fit, type = "score", weighted=F))
    score <- score[subco.noncase,,drop=F] ## Scores for controls
    
    nn <- as.vector(cohort.size)
    d <- table(factor(stratum[cases], levels = sort(unique(stratum))))
    nn0 <- nn - as.vector(d)
    m0 <- table(stratum[subco.noncase])
    wgt <- as.vector(nn0/m0)
    
    st <- stratum[subco.noncase] ## Stratum indicators for controls
    kk    <- ncol(score)
    sto <- st %o% rep(1,kk)
    Index <- col(score)
    tscore <- tapply(score,list(sto,Index),mean) ## Within stratum control score means
    pscore <- tapply(score,list(sto,Index))
    score <- score-tscore[pscore] ## Subtract off within stratum score means
    delta <- matrix(0,kk,kk)
    opt <- NULL
    for (j in 1:nstrata) {
      temp <- t(score[st==j,])%*%score[st==j,]/(m0[j]-1) ## Borgan equation (19)
      delta <- delta+(wgt[j]-1)*nn0[j]*temp ## Borgan equation (17)
      if(is.null(opt)) 
        opt <- nn0[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var)) 
      else
        opt <- rbind(opt, nn0[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var))) 
    }
    z <- apply(opt,2,sum)
    fit$opt <- sweep(opt,2,z,FUN="/")
    fit$phase2var<-fit$naive.var %*% delta %*% fit$naive.var
    fit$naive.var <- fit$naive.var+fit$phase2var
    fit$var<-fit$naive.var
    return(fit)
  }
  else if (cch.method == "Samuelson2007"){
    wgt <- df$weights
    subco.cases <- split(subco.cases, df[[stratum]][subco.cases])
    cases <- split(cases, df[[stratum]][cases])
    subco.noncase <- split(subco.noncase, df[[stratum]][subco.noncase])
    
    N_str <- cohort.size
    D_str <- lengths(cases)
    N_nc <- N_str - D_str
    mn1 <- lengths(subco.noncase)
    db <- as.matrix(resid(fit, type = "dfbeta"))
    
    phase2_terms <- lapply(seq_along(N_str), function(j) {
      idx_j <- subco.noncase[[j]]
      db_j  <- db[idx_j, , drop = FALSE]
      dbm   <- colMeans(db_j)
      db0   <- sweep(db_j, 2, dbm)
      
      bessel <- mn1[j] / (mn1[j] - 1)
      fpc <- 1 - mn1[j] / N_nc[j]
      
      bessel * fpc * crossprod(db0)
    })
    
    fit$phase2var <- Reduce(`+`, phase2_terms)
    fit$naive.var <- fit$naive.var + fit$phase2var
    fit$var <- fit$naive.var
    return(fit)
  }
}