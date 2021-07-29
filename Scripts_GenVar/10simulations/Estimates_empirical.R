
#######################################################################
#         Genetic and genic variances obtained from AlphaBayes        #
#             Letícia Lara, Ivan Pocrnic, Gregor Gorjanc              #
#                             June, 2020                              #
#######################################################################


library(AlphaSimR)
library(dplyr)


#### Loading and manipulating the data ####
load(paste0("../AlphaSimR/Results/Burnin", REP, ".RData"))
load(paste0("../AlphaSimR/Results/all_stages_and_years", REP, ".RData"))

emp.sample <- matrix(as.numeric(scan(paste0("Sample", REP, "_empirical.txt"), what = "numeric")), nrow = 10500)
emp.sample <- t(emp.sample)
dim(emp.sample)

sample.effect <- lapply(seq_len(nrow(emp.sample)), function(i) emp.sample[i,])
length(sample.effect)


#### Predictions of each stage and year of the breeding program ####
vector.EYT <- rep(c(rep(1, 10), rep(2, 10)), 6)
EYT1 <- EYTc[vector.EYT == 1]
EYT2 <- EYTc[vector.EYT == 2]

for(year in 16:21){
  load(paste0("../AlphaSimR/Results/stages_year", year, REP, ".RData"))

  # Parents - 70 x 10500
  vector.year <- rep(16:21, each=70)
  Parents.year <- parents.c[vector.year == year]
  Parents.Genetic <- lapply(sample.effect, function(x) stage.year.Parents$snp %*% x)  # Wα
  Parents.Genetic.est <- sapply(Parents.Genetic, function(x) var(x))
  Parents.Genetic.var <- Parents.Genetic.est*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)
  Parents.pq <- apply(stage.year.Parents$snp, 2, function(x) var(x))     
  Parents.Genic <- sapply(sample.effect, function(x) Parents.pq %*% x^2)   # 2pqα^2
  Parents.Genic.var <- Parents.Genic*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)

  # F1 - 100 x 10500
  vector.year <- rep(16:21, each=100)
  F1.year <- F1c[vector.year == year]
  F1.Genetic <- lapply(sample.effect, function(x) stage.year.F1$snp %*% x)  # Wα 
  F1.Genetic.est <- sapply(F1.Genetic, function(x) var(x))
  F1.Genetic.var <- F1.Genetic.est*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)
  F1.pq <- apply(stage.year.F1$snp, 2, function(x) var(x))     
  F1.Genic <-  sapply(sample.effect, function(x) F1.pq %*% x^2)   # 2pqα^2
  F1.Genic.var <- F1.Genic*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)

  # HDRW - 10000 x 10500 from DH and 500 x 10500 from GS
  vector.year <- c(rep(16:20, each=10000), rep(21, 500))
  HDRW.year <- HDRWc[vector.year == year]
  HDRW.Genetic <- lapply(sample.effect, function(x) stage.year.HDRW$snp %*% x)  # Wα
  HDRW.Genetic.est <- sapply(HDRW.Genetic, function(x) var(x))
  HDRW.Genetic.var <- HDRW.Genetic.est*(nrow(stage.year.HDRW$snp)-1)/nrow(stage.year.HDRW$snp)
  HDRW.pq <- apply(stage.year.HDRW$snp, 2, function(x) var(x))     
  HDRW.Genic <- sapply(sample.effect, function(x) HDRW.pq %*% x^2)   # 2pqα^2
  HDRW.Genic.var <- HDRW.Genic*(nrow(stage.year.HDRW$snp)-1)/nrow(stage.year.HDRW$snp)

  # PYT - 500 x 10500
  vector.year <- rep(16:21, each=500)
  PYT.year <- PYTc[vector.year == year]
  PYT.Genetic <- lapply(sample.effect, function(x) stage.year.PYT$snp %*% x)  # Wα
  PYT.Genetic.est <- sapply(PYT.Genetic, function(x) var(x))
  PYT.Genetic.var <- PYT.Genetic.est*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)
  PYT.pq <- apply(stage.year.PYT$snp, 2, function(x) var(x))     
  PYT.Genic <- sapply(sample.effect, function(x) PYT.pq %*% x^2)   # 2pqα^2
  PYT.Genic.var <- PYT.Genic*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)

  # AYT - 50 x 10500
  vector.year <- rep(16:21, each=50)
  AYT.year <- AYTc[vector.year == year]
  AYT.Genetic <- lapply(sample.effect, function(x) stage.year.AYT$snp %*% x)  # Wα
  AYT.Genetic.est <- sapply(AYT.Genetic, function(x) var(x))
  AYT.Genetic.var <- AYT.Genetic.est*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)
  AYT.pq <- apply(stage.year.AYT$snp, 2, function(x) var(x))     
  AYT.Genic <- sapply(sample.effect, function(x) AYT.pq %*% x^2)   # 2pqα^2
  AYT.Genic.var <- AYT.Genic*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)

  # EYT - 20 x 10500
  vector.year <- rep(16:21, each=20)
  EYT.year <- EYTc[vector.year == year]
  EYT.Genetic <- lapply(sample.effect, function(x) stage.year.EYT$snp %*% x)  # Wα
  EYT.Genetic.est <- sapply(EYT.Genetic, function(x) var(x))
  EYT.Genetic.var <- EYT.Genetic.est*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)
  EYT.pq <- apply(stage.year.EYT$snp, 2, function(x) var(x))     
  EYT.Genic <- sapply(sample.effect, function(x) EYT.pq %*% x^2)   # 2pqα^2
  EYT.Genic.var <- EYT.Genic*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)

  # EYT1 - 10 x 10500
  snp.EYT1 <- stage.year.EYT$snp[1:10,]
  vector.year <- rep(16:21, each=10)
  EYT1.year <- EYT1[vector.year == year]
  EYT1.true.Genetic <- varG(EYT1.year)
  EYT1.true.Genic <- genicVarG(EYT1.year)
  EYT1.Genetic <- lapply(sample.effect, function(x) stage.year.EYT1$snp %*% x)  # Wα
  EYT1.Genetic.est <- sapply(EYT1.Genetic, function(x) var(x))
  EYT1.Genetic.var <- EYT1.Genetic.est*(nrow(stage.year.EYT1$snp)-1)/nrow(stage.year.EYT1$snp)
  EYT1.pq <- apply(stage.year.EYT1$snp, 2, function(x) var(x))     
  EYT1.Genic <- sapply(sample.effect, function(x) EYT1.pq %*% x^2)   # 2pqα^2
  EYT1.Genic.var <- EYT1.Genic*(nrow(stage.year.EYT1$snp)-1)/nrow(stage.year.EYT1$snp)

  # EYT2 - 10 x 10500
  snp.EYT2 <- stage.year.EYT$snp[11:20,]
  vector.year <- rep(16:21, each=10)
  EYT2.year <- EYT2[vector.year == year]
  EYT2.true.Genetic <- varG(EYT2.year)
  EYT2.true.Genic <- genicVarG(EYT2.year)
  EYT2.Genetic <- lapply(sample.effect, function(x) stage.year.EYT2$snp %*% x)  # Wα
  EYT2.Genetic.est <- sapply(EYT2.Genetic, function(x) var(x))
  EYT2.Genetic.var <- EYT2.Genetic.est*(nrow(stage.year.EYT2$snp)-1)/nrow(stage.year.EYT2$snp)
  EYT2.pq <- apply(stage.year.EYT2$snp, 2, function(x) var(x))     
  EYT2.Genic <- sapply(sample.effect, function(x) EYT2.pq %*% x^2)   # 2pqα^2
  EYT2.Genic.var <- EYT2.Genic*(nrow(stage.year.EYT2$snp)-1)/nrow(stage.year.EYT2$snp)

  save(Parents.Genetic.var, Parents.Genic.var, F1.Genetic.var, F1.Genic.var, 
       HDRW.Genetic.var, HDRW.Genic.var, PYT.Genetic.var, PYT.Genic.var, 
       AYT.Genetic.var, AYT.Genic.var, EYT.Genetic.var, EYT.Genic.var, 
       EYT1.true.Genetic, EYT1.true.Genic, EYT1.Genetic.var, EYT1.Genic.var, 
       EYT2.true.Genetic, EYT2.true.Genic, EYT2.Genetic.var, EYT2.Genic.var, 
       file = paste0("./Simulation", simulation, "/empiricalBayes/estimates.empiricalBayes_year", year, REP, ".RData"))
}






