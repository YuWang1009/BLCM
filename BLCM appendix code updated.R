# first install jags before running model, available at https://mcmc-jags.sourceforge.io/

rm(list = ls())
dev.off()
gc()

getwd()
setwd("...")

# install.packages("R2jags")
library(R2jags)

model.text1 <- "model{
for (i in 1:651){
  logit(p[i]) <- alpha[region[i]] + beta_age[age[i]] + U[farm[i]]               # predict probability to be seropositive for cattle i
  z[i] ~ dbern(p[i])                                                            # sample the seropositivity status for cattle i
  pRBT[i] <- z[i]*seRBT+(1-z[i])*(1-spRBT)                                      # apparent RBT prevalence
  RBT[i] ~ dbern(pRBT[i])                                                       # sample the apparent RBT result
  TRBT[i] <- RBT[i] + 1 # index for RBT status                                  # index for RBT status, TRBT = 2 RBT+; TRBT = 1, RBT-
  pSAT[i] <- z[i]*seSATc[TRBT[i]]+(1-z[i])*(1-spSATc[TRBT[i]])                  # apparent SAT prevalence expressed conditionally on RBT
  SAT[i] ~ dbern(pSAT[i]) # index for SAT status                                # sample the apparent SAT result
  dummySATRBTZ[i] <- 1+z[i]+2*RBT[i]+4*SAT[i]                                   # index for eight combination by disease status, RBT, SAT
  # 1 if true status=0=RBT=SAT
  # 2 if true status=1, and 0=RBT=SAT (both TN)
  # 3 if true status=0, RBT=1 (FP), SAT=0 (TN)
  # 4 if true status=1, RBT=1 (TP), SAT=0 (FN)
  # 5 if true status=0, RBT=0 (TN), SAT=1 (FP)
  # 6 if true status=1, RBT=0 (FN), SAT=1 (TP)
  # 7 if true status=0, RBT=1 (FP), SAT=1 (FP)
  # 8 if true status=1, RBT=1 (TP), SAT=1 (TP)

  mP[i] ~ dnorm(mu.mP[dummySATRBTZ[i]], prec.mP[dummySATRBTZ[i]])               # Gassusian mixture model for FPA (mP), log-transformed

  PI[i] ~ dnorm(mu.cor.PI[i, dummySATRBTZ[i]], prec.mPPI[dummySATRBTZ[i]])      # Gassusian mixture model for ELISA (PI), log-transformed
  for(j in 1:8){
   # this is an artefact to create correlation. Has no interpretation for us
   mu.cor.PI[i, j] <- mu.PI[j]+rho[j]*sqrt(var.PI[j]/var.mP[j])*(mP[i]-mu.mP[j])
  }
}

  # proportion of each test status
  pp[1] <- (1-seRBT)*(1-seSATc[1])  # true status = 1, RBT = 0, SAT= 0; 0 - negative, 1 - positive
  pp[2] <- seRBT*(1-seSATc[2])      # true status = 1, RBT = 1, SAT= 0
  pp[3] <- (1-seRBT)*seSATc[1]      # true status = 1, RBT = 0, SAT= 1
  pp[4] <- seRBT*seSATc[2]          # true status = 1, RBT = 1, SAT= 1

  pn[1] <- spRBT*spSATc[1]          # true status = 0, RBT = 0, SAT= 0
  pn[2] <- (1-spRBT)*spSATc[2]      # true status = 0, RBT = 1, SAT= 0
  pn[3] <- spRBT*(1-spSATc[1])      # true status = 0, RBT = 0, SAT= 1
  pn[4] <- (1-spRBT)*(1-spSATc[2])  # true status = 0, RBT = 1, SAT= 1

  # derived parameters
  seSAT <- seSATc[1]*(1-seRBT)+seSATc[2]*seRBT
  spSAT <- spSATc[1]*spRBT+spSATc[2]*(1-spRBT)
  rhoDRBTSAT <- seRBT*(seSATc[2]-seSAT)/sqrt(seRBT*seSAT*(1-seRBT)*(1-seSAT))     # rho between RBT and SAT given D+
  rhoNRBTSAT <- (1-spRBT)*(spSAT-spSATc[2])/sqrt(spRBT*spSAT*(1-spRBT)*(1-spSAT)) # rho between RBT and SAT given D-

  muDFPA <- mu.mP[2]*pp[1]+mu.mP[4]*pp[2]+mu.mP[6]*pp[3]+mu.mP[8]*pp[4]     # mixture mean FPA value for diseased animals (D+)
  muDRBTFPA <- mu.mP[4]*pp[2]+mu.mP[8]*pp[4]                                # mixture mean FPA value given RBT+, D+
  muDSATFPA <- mu.mP[6]*pp[3]+mu.mP[8]*pp[4]                                # mixture mean FPA value given SAT+, D+
  varDFPA <- mu.mP[2]*mu.mP[2]*pp[1]+mu.mP[4]*mu.mP[4]*pp[2]+mu.mP[6]*mu.mP[6]*pp[3]+mu.mP[8]*mu.mP[8]*pp[4]+
  var.mP[2]*pp[1]+var.mP[4]*pp[2]+var.mP[6]*pp[3]+var.mP[8]*pp[4]-muDFPA*muDFPA      # variance of FPA for diseased animals
  rhoDRBTFPA <- (muDRBTFPA-muDFPA*seRBT)/sqrt(seRBT*(1-seRBT)*varDFPA)      # rho between RBT and FPA given D+
  rhoDSATFPA <- (muDSATFPA-muDFPA*seSAT)/sqrt(seSAT*(1-seSAT)*varDFPA)      # rho between SAT adn FPA given D+

  muNFPA <- mu.mP[1]*pn[1]+mu.mP[3]*pn[2]+mu.mP[5]*pn[3]+mu.mP[7]*pn[4]     # mixture mean FPA value for healthy animals (D-)
  muNRBTFPA <- mu.mP[3]*pn[2]+mu.mP[7]*pn[4]                                # mixture mean FPA value given RBT-, D-
  muNSATFPA <- mu.mP[5]*pn[3]+mu.mP[7]*pn[4]                                # mixture mean FPA value given SAT-, D+
  varNFPA <- mu.mP[1]*mu.mP[1]*pn[1]+mu.mP[3]*mu.mP[3]*pn[2]+mu.mP[5]*mu.mP[5]*pn[3]+mu.mP[7]*mu.mP[7]*pn[4]+
  var.mP[1]*pn[1]+var.mP[3]*pn[2]+var.mP[5]*pn[3]+var.mP[7]*pn[4]-muNFPA*muNFPA       # variance of FPA for D- animals
  rhoNRBTFPA <- (muNRBTFPA-muNFPA*(1-spRBT))/sqrt(spRBT*(1-spRBT)*varNFPA)  # rho between RBT and FPA given D-
  rhoNSATFPA <- (muNSATFPA-muNFPA*(1-spSAT))/sqrt(spSAT*(1-spSAT)*varNFPA)  # rho between SAT adn FPA given D-

  # same as FPA
  muDELISA <- mu.PI[2]*pp[1]+mu.PI[4]*pp[2]+mu.PI[6]*pp[3]+mu.PI[8]*pp[4]
  muDRBTELISA <- mu.PI[4]*pp[2]+mu.PI[8]*pp[4]
  muDSATELISA <- mu.PI[6]*pp[3]+mu.PI[8]*pp[4]
  varDELISA <- mu.PI[2]*mu.PI[2]*pp[1]+mu.PI[4]*mu.PI[4]*pp[2]+mu.PI[6]*mu.PI[6]*pp[3]+mu.PI[8]*mu.PI[8]*pp[4]+
  var.PI[2]*pp[1]+var.PI[4]*pp[2]+var.PI[6]*pp[3]+var.PI[8]*pp[4]-muDELISA*muDELISA
  rhoDRBTELISA <- (muDRBTELISA-muDELISA*seRBT)/sqrt(seRBT*(1-seRBT)*varDELISA)
  rhoDSATELISA <- (muDSATELISA-muDELISA*seSAT)/sqrt(seSAT*(1-seSAT)*varDELISA)

  muNELISA <- mu.PI[1]*pn[1]+mu.PI[3]*pn[2]+mu.PI[5]*pn[3]+mu.PI[7]*pn[4]
  muNRBTELISA <- mu.PI[3]*pn[2]+mu.PI[7]*pn[4]
  muNSATELISA <- mu.PI[5]*pn[3]+mu.PI[7]*pn[4]
  varNELISA <- mu.PI[1]*mu.PI[1]*pn[1]+mu.PI[3]*mu.PI[3]*pn[2]+mu.PI[5]*mu.PI[5]*pn[3]+mu.PI[7]*mu.PI[7]*pn[4]+
  var.PI[1]*pn[1]+var.PI[3]*pn[2]+var.PI[5]*pn[3]+var.PI[7]*pn[4] - muNELISA*muNELISA
  rhoNRBTELISA <- (muNRBTELISA-muNELISA*(1-spRBT))/sqrt(spRBT*(1-spRBT)*varNELISA)
  rhoNSATELISA <- (muNSATELISA-muNELISA*(1-spSAT))/sqrt(spSAT*(1-spSAT)*varNELISA)

  for (k in 1:100){
  c[k] <- -0.3+k*0.005  # cut-off range (74mP - 122mP)
  seFPA[k] <- 1-(pp[1]*phi((c[k]-mu.mP[2])/sqrt(var.mP[2]))+pp[2]*phi((c[k]-mu.mP[4])/sqrt(var.mP[4]))+
  pp[3]*phi((c[k]-mu.mP[6])/sqrt(var.mP[6]))+pp[4]*phi((c[k]-mu.mP[8])/sqrt(var.mP[8])))
  fpFPA[k] <- 1-(pn[1]*phi((c[k]-mu.mP[1])/sqrt(var.mP[1]))+pn[2]*phi((c[k]-mu.mP[3])/sqrt(var.mP[3]))+
  pn[3]*phi((c[k]-mu.mP[5])/sqrt(var.mP[5]))+pn[4]*phi((c[k]-mu.mP[7])/sqrt(var.mP[7])))
  spFPA[k] <- 1-fpFPA[k]
  youdenFPA[k] <- seFPA[k]+spFPA[k]

  d[k] <- -1.5+k*0.01   # cut-off range (0.22,0.60)
  seELISA[k] <- 1-(pp[1]*phi((d[k]-mu.PI[2])/sqrt(var.PI[2]))+pp[2]*phi((d[k]-mu.PI[4])/sqrt(var.PI[4]))+
  pp[3]*phi((d[k]-mu.PI[6])/sqrt(var.PI[6]))+pp[4]*phi((d[k]-mu.PI[8])/sqrt(var.PI[8])))
  fpELISA[k] <- 1-(pn[1]*phi((d[k]-mu.PI[1])/sqrt(var.PI[1]))+pn[2]*phi((d[k]-mu.PI[3])/sqrt(var.PI[3]))+
  pn[3]*phi((d[k]-mu.PI[5])/sqrt(var.PI[5]))+pn[4]*phi((d[k]-mu.PI[7])/sqrt(var.PI[7])))
  spELISA[k] <- 1-fpELISA[k]
  youdenELISA[k] <- seELISA[k]+spELISA[k]
  }

  # Parallel and series interpretation
  se.RBT.SAT.ser <- pp[4]
  sp.RBT.SAT.ser <- 1-pn[4]
  se.RBT.SAT.par <- 1-pp[1]
  sp.RBT.SAT.par <- pn[1]

  # priors
  p1 ~ dbeta(1.14, 14.56) # 95% certain prevalence less than 0.2 at 0.01 mode for region 1
  p2 ~ dbeta(1.06,1.32)   # 95% certain prevalence great than 0.01 at 0.15 mode for region 2

  # logit-transforming prior distributions of region-specific prevalence
  alpha[1] <- logit(p1)
  alpha[2] <- logit(p2)

  beta_age[1] <- 1
  for (i in 2:5){
   beta_age[i] ~ dnorm(0, 0.1)
  }
  for (i in 1:6){
   U[i] ~ dnorm(0, ivar)
  }
  ivar <- 1/sd/sd
  sd ~ dunif(0, 1)

  # sensitivity and specificity for RBT
  seRBT ~ dbeta(13.58,3.91)T(1-spRBT, ) # 95% certain se more than 0.6 at 0.812 mode
  spRBT ~ dbeta(9.92,2.42)              # 95% certain sp more than 0.6 at 0.863 mode

  # conditional sensitivity and specificity for SAT on RBT
  seSATc[1] ~ dbeta(21.40,7.48)T(1-spSATc[1], )   # 95% certain se more than 0.6 at 0.759 mode
  seSATc[2] ~ dbeta(21.40,7.48)T(1-spSATc[2], )
  spSATc[1] ~ dbeta(9.92,2.42)                    # 95% certain sp more than 0.3 at 0.863 mode
  spSATc[2] ~ dbeta(9.92,2.42)

  # priors for mP, 8 combinations, set some limitations as diseased animals generally have higher antibody titers than non-diseased animals
  for (i in 1:4){
    mu.mP[2*i-1] ~ dunif(-0.3,0.0)         # mean of non-diseased status (D-), raw negative mP generally ranges from 70 to 100, so the log(mP) ranges from 4 to 5
    mu.mP[2*i] ~ dunif(mu.mP[2*i-1],0.6)   # mean of diseased status (D+)
    mu.PI[2*i-1] ~ dunif(-2.3,-0.7)
    mu.PI[2*i] ~ dunif(mu.PI[2*i-1],0.3)
  }

  for (i in 1:8){
    prec.mP[i] ~ dgamma(0.01,0.01)
    var.mP[i] <- 1/prec.mP[i]
    prec.PI[i] ~ dgamma(0.01,0.01)
    var.PI[i] <- 1/prec.PI[i]
    rho[i] ~ dunif(0,1) # force to positive correlation
    prec.mPPI[i] <- prec.PI[i]/(1-rho[i]^2)
  }
}"

data <- read.csv("test data.csv")
dat.test <- list(age = data$parity+1, farm = data$farm, region = data$region,
                 RBT = data$RBT, SAT = data$SAT, mP = log(data$mP), PI = log(data$PI))


params <- c("alpha","beta_age","U","z","seRBT", "spRBT","seSAT", "spSAT",
            "seSATc", "spSATc", "rho","seFPA","spFPA","rhoDRBTFPA","rhoNRBTFPA","rhoDSATFPA","rhoNSATFPA",
            "rhoDRBTELISA","rhoNRBTELISA","rhoDSATELISA","rhoNSATELISA",
            "youdenFPA","seELISA","spELISA","youdenELISA","mu.mP","mu.PI","prec.mP","prec.PI",
            "se.RBT.SAT.ser","sp.RBT.SAT.ser","se.RBT.SAT.par","sp.RBT.SAT.par")

start <- Sys.time()
fit <- jags(
  data = dat.test,
  parameters.to.save = params,
  model.file = textConnection(model.text1),
  n.chains = 2,
  n.iter = 20000,
  n.burnin = 5000,
  n.thin = 1,
  DIC = TRUE
)

end <- Sys.time()
(end-start)  # it will run about 2 hours

fit

traceplot(fit)

results <- as.data.frame(fit[["BUGSoutput"]][["summary"]])
results

