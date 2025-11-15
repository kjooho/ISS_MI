library(mice)
library(survival)
library(mitools)
source("R/supersampling.R")
source("R/addph2var.R")
source("R/smcfcs.cc.R")
source("R/generate.data.R")

ccdesign <- "case.cohort"  # "case.cohort" or "stratified.case.cohort"
dat_list <- generate.data(0.99, design = ccdesign, interaction = F)
data_mi <- dat_list$dat.mi
data_smc <- dat_list$dat.smc
CC_sample <- data_smc[data_smc$ind.ph2==1, ]

smformula = "Surv(eventime, ind.fail) ~ z1 + z2 + z3 + xc + xb"
sub.smformula = "Surv(eventime, ind.fail) ~ z1 + z2 + z3"
in.subco <- "ind.ini"
stratum <- "stratum"

if (ccdesign == "case.cohort") {
  cohort.size <- 25000
  n.super <- 750
  cch.method  <- "LinYing"
} else if (ccdesign == "stratified.case.cohort") {
  cohort.size <- n.super <- table(data_mi[["stratum"]])
  n.super[1:6] <- 125
  cch.method  <- "Samuelson2007"
}

# Method and predictor matrix
meth_mice <- make.method(data_mi)
meth_mice["xc"] <- "norm"
meth_mice["xb"] <- "logreg"

meth_smc <- make.method(data_smc)
meth_smc["xc"] <- "norm"
meth_smc["xb"] <- "logreg"

# remember to exclude weights in the imputation model!
pred_mice <- make.predictorMatrix(data_mi)
if (ccdesign == "case.cohort"){
  pred_mice[c("xc","xb"), c("eventime","weights","id","ind.ph2")] <- 0
} else if (ccdesign == "stratified.case.cohort"){
  pred_mice[c("xc","xb"), c("eventime","weights","id","ind.ph2","stratum")] <- 0
}

pred_smc <- make.predictorMatrix(data_smc)
if (ccdesign == "case.cohort"){
  pred_smc[c("xc","xb"), c("eventime","weights","id","ind.ini","ind.fail","ind.ph2")] <- 0
} else if (ccdesign == "stratified.case.cohort"){
  pred_smc[c("xc","xb"), c("eventime","weights","id","ind.ini","ind.fail","ind.ph2","stratum")] <- 0
}

##################################################################
#                              MICE                              #
##################################################################
# Divide observations by cases and censored units
ind.fail <- as.character(as.formula(smformula)[[2]][[3]])
subc_noncase <- data_mi[data_mi[,in.subco] == 1 | data_mi[,ind.fail] == 0, ]
outsubc_case <- data_mi[data_mi[,in.subco] == 0 & data_mi[,ind.fail] == 1, ]
outsubc_case_dt <- data.table::as.data.table(outsubc_case)

# Impute using censored units
complete_list <- subc_noncase %>% 
  mice(m=10, numit = 50, pred=pred_mice, method=meth_mice, printFlag=F) %>%
  mice::complete("all", include = F)

# Add cases to imputed data and conduct case-cohort analysis
fits_mi <- lapply(complete_list, function(df) {
  df_dt <- data.table::as.data.table(df)
  full_dt <- data.table::rbindlist(list(df_dt, outsubc_case_dt), use.names = TRUE, fill = TRUE)
  survival::coxph(as.formula(smformula), control = survival::coxph.control(timefix = FALSE), data = full_dt, cluster = full_dt$id)
})

summary(pool(fits_mi))

##################################################################
#                            MICE ISS                            #
##################################################################
super <- supersampling(data=data_mi, sub.smformula=sub.smformula, in.subco=in.subco, 
                        n.super=n.super, random=F, design="case.cohort", stratum=stratum)
supersample <- data_mi[super$indices, , drop = FALSE]
supersample[["weights"]] <- super$weights

ind.fail <- as.character(as.formula(smformula)[[2]][[3]])
subc_noncase <- supersample[supersample[,in.subco] == 1 | supersample[,ind.fail] == 0, ]
outsubc_case <- supersample[supersample[,in.subco] == 0 & supersample[,ind.fail] == 1, ]
outsubc_case_dt <- data.table::as.data.table(outsubc_case)

# Impute using subcohort units
complete_1 <- subc_noncase %>%
  mice(m=10, numit = 50, pred=pred_mice, method=meth_mice, printFlag=F) %>%
  mice::complete("all", include = F)

cohort.size <- nrow(data_mi)
stratified.cohort.size <- table(data_mi[[stratum]])

# Add non-subcohort cases to imputed data and conduct case-cohort analysis
fits_miiss <- lapply(complete_1, function(df) {
  df_dt <- data.table::as.data.table(df)
  full_dt <- data.table::rbindlist(list(df_dt, outsubc_case_dt), use.names = TRUE, fill = TRUE)
  idx_nc <- full_dt[[ind.fail]] == 0L
  full_dt[idx_nc, (in.subco) := 1L]
  fit <- survival::coxph(as.formula(smformula), 
                         data = full_dt, 
                         control = survival::coxph.control(timefix = FALSE),
                         cluster = full_dt$id,
                         weights = full_dt$weights)
  if (cch.method=="LinYing"){
    fit <- add.ph2var(fit=fit, df=full_dt, smformula=smformula, in.subco=in.subco, cohort.size=cohort.size, cch.method=cch.method, robust=F)
  } else if (cch.method=="Samuelson2007"){
    fit <- add.ph2var(fit=fit, df=full_dt, smformula=smformula, in.subco=in.subco, cohort.size=stratified.cohort.size, stratum=stratum, cch.method=cch.method)
  }
  fit
})
summary(MIcombine(fits_miiss))

##################################################################
#                              SMCFCS                            #
##################################################################
imp.rj <- smcfcs.cc(data_smc, smtype="coxph", smformula=smformula, pred=pred_smc, method=meth_smc,
                    m=10, numit=50, rjlimit=10000, in.subco=in.subco, super=F)
imp.obj <- imputationList(imp.rj$impDatasets)
fits_smc <- with(imp.obj, survival::coxph(as.formula(smformula), control = survival::coxph.control(timefix = FALSE), cluster=id))
summary(MIcombine(fits_smc))

##################################################################
#                            SMCFCS ISS                          #
##################################################################
super <- supersampling(data=data_smc, sub.smformula=sub.smformula, in.subco=in.subco, 
                        n.super=n.super, random=F, design="case.cohort", stratum = stratum)
supersample <- data_smc[super$indices, , drop = FALSE]
supersample[["weights"]] <- super$weights

imp.rj <- smcfcs.cc(supersample, smtype="coxph", smformula=smformula, pred=pred_smc, 
                    method=meth_smc, m=10, numit=50, rjlimit=10000, in.subco=in.subco, super=T)
imp_dat <- imp.rj$impDatasets

cohort.size <- nrow(data_smc)
stratified.cohort.size <- table(data_smc[[stratum]])
ind.fail <- as.character(as.formula(smformula)[[2]][[3]])

fits_smciss <- lapply(imp_dat, function(df) {
  fit <- survival::coxph(as.formula(smformula), 
                         data = df, 
                         control = survival::coxph.control(timefix = FALSE),
                         cluster = df$id,
                         weights = df$weights)
  df[df[[ind.fail]]==0, in.subco] <- 1
  if (cch.method=="LinYing"){
    fit <- add.ph2var(fit=fit, df=df, smformula=smformula, in.subco=in.subco, cohort.size=cohort.size, cch.method=cch.method, robust=F)
  } else if (cch.method=="Samuelson2007"){
    fit <- add.ph2var(fit=fit, df=df, smformula=smformula, in.subco=in.subco, cohort.size=stratified.cohort.size, stratum=stratum, cch.method=cch.method)
  }
  fit
})
summary(MIcombine(fits_smciss))

