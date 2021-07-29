
#######################################################################
#         Genetic and genic variances obtained from AlphaBayes        #
#                         Using SVD approach                          #
#             Letícia Lara, Ivan Pocrnic, Gregor Gorjanc              #
#                             March, 2019                             #
#######################################################################

#options(repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))
#install.packages("INLA", dep = TRUE)

library(AlphaSimR)
library(tidyverse)
library(irlba)
library(dplyr)

load("../../Burnin.RData")
load("../../all_stages_and_year.RData")
load("./prcomp.irlba_10pc.RData")

eigenvectors <- irlba_10pc$rotation[,1:10]

sample.est <-  matrix(as.numeric(scan("./MarkerSamples_10pc.txt", what = "numeric")), nrow = 10)
pc.mat <- t(sample.est)

set.seed(24082019)
sample <- sort(sample(1:90000, 1000))
pc.sample <- pc.mat[sample,]
dim(pc.sample)

sample.effect <- lapply(seq_len(nrow(pc.sample)), function(i) pc.sample[i,])
str(sample.effect[1:3])

save(sample.effect, file = "sample.effect_10pc.RData")

vector.EYT <- rep(c(rep(1,10), rep(2,10)),6)
EYT1 <- EYTc[vector.EYT == 1]
EYT2 <- EYTc[vector.EYT == 2]

for(year in 16:21){
  load(paste0("../../output_sim/stages_year",year,".RData"))
  
  # Parents - 70 x 10
  vector.year <- rep(16:21, each=70)
  Parents.year <- parents.c[vector.year == year]
  Parents.Genetic <- lapply(sample.effect, function(x) stage.year.Parents$snp %*% eigenvectors %*% x)  # Wvα*
  Parents.Genetic.est <- sapply(Parents.Genetic, function(x) var(x))
  Parents.Genetic.var <- Parents.Genetic.est*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)
  Parents.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  Parents.pq <- apply(stage.year.Parents$snp, 2, function(x) var(x))     
  Parents.Genic <- sapply(Parents.alpha, function(x) Parents.pq %*% x^2)   # 2pqα^2
  Parents.Genic.var <- Parents.Genic*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)
  Parents.randGenic.true <- genParam(Parents.year)$randGenicVarG
  
  # F1 - 100 x 10
  vector.year <- rep(16:21, each=100)
  F1.year <- F1c[vector.year == year]
  F1.Genetic <- lapply(sample.effect, function(x) stage.year.F1$snp %*% eigenvectors %*% x)  # Wvα*
  F1.Genetic.est <- sapply(F1.Genetic, function(x) var(x))
  F1.Genetic.var <- F1.Genetic.est*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)
  F1.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  F1.pq <- apply(stage.year.F1$snp, 2, function(x) var(x))     
  F1.Genic <- sapply(F1.alpha, function(x) F1.pq %*% x^2)   # 2pqα^2
  F1.Genic.var <- F1.Genic*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)
  F1.randGenic.true <- genParam(F1.year)$randGenicVarG

  # HDRW - 10000 x 10 from DH and 500 x 10 from GS
  if(year != 21){
    seed_value <- as.integer(paste0("2019",year))*2
    set.seed(seed_value)
    sampled <- sample(1:10000, 1000, replace = F)
    sampled.HDRW.snp <- stage.year.HDRW$snp[sort(sampled),]  # 1000 x 10500
  }else{
    sampled.HDRW.snp <- stage.year.HDRW$snp  # 500 x 10500
  }
  vector.year <- c(rep(16:20, each=10000), rep(21, 500))
  HDRW.year <- HDRWc[vector.year == year]
  HDRW.Genetic <- lapply(sample.effect, function(x) sampled.HDRW.snp %*% eigenvectors %*% x)  # Wvα*
  HDRW.Genetic.est <- sapply(HDRW.Genetic, function(x) var(x))
  HDRW.Genetic.var <- HDRW.Genetic.est*(nrow(sampled.HDRW.snp)-1)/nrow(sampled.HDRW.snp)
  HDRW.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  HDRW.pq <- apply(sampled.HDRW.snp, 2, function(x) var(x))     
  HDRW.Genic <- sapply(HDRW.alpha, function(x) HDRW.pq %*% x^2)   # 2pqα^2
  HDRW.Genic.var <- HDRW.Genic*(nrow(sampled.HDRW.snp)-1)/nrow(sampled.HDRW.snp)
  HDRW.randGenic.true <- genParam(HDRW.year)$randGenicVarG
  
  # PYT - 500 x 10
  vector.year <- rep(16:21, each=500)
  PYT.year <- PYTc[vector.year == year]
  PYT.Genetic <- lapply(sample.effect, function(x) stage.year.PYT$snp %*% eigenvectors %*% x)  # Wvα*
  PYT.Genetic.est <- sapply(PYT.Genetic, function(x) var(x))
  PYT.Genetic.var <- PYT.Genetic.est*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)
  PYT.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  PYT.pq <- apply(stage.year.PYT$snp, 2, function(x) var(x))     
  PYT.Genic <- sapply(PYT.alpha, function(x) PYT.pq %*% x^2)   # 2pqα^2
  PYT.Genic.var <- PYT.Genic*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)
  PYT.randGenic.true <- genParam(PYT.year)$randGenicVarG
  
  # AYT - 50 x 10
  vector.year <- rep(16:21, each=50)
  AYT.year <- AYTc[vector.year == year]
  AYT.Genetic <- lapply(sample.effect, function(x) stage.year.AYT$snp %*% eigenvectors %*% x)  # Wvα*
  AYT.Genetic.est <- sapply(AYT.Genetic, function(x) var(x))
  AYT.Genetic.var <- AYT.Genetic.est*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)
  AYT.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  AYT.pq <- apply(stage.year.AYT$snp, 2, function(x) var(x))     
  AYT.Genic <- sapply(AYT.alpha, function(x) AYT.pq %*% x^2)   # 2pqα^2
  AYT.Genic.var <- AYT.Genic*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)
  AYT.randGenic.true <- genParam(AYT.year)$randGenicVarG
  
  # EYT - 20 x 10
  vector.year <- rep(16:21, each=20)
  EYT.year <- EYTc[vector.year == year]
  EYT.Genetic <- lapply(sample.effect, function(x) stage.year.EYT$snp %*% eigenvectors %*% x)  # Wvα*
  EYT.Genetic.est <- sapply(EYT.Genetic, function(x) var(x))
  EYT.Genetic.var <- EYT.Genetic.est*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)
  EYT.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  EYT.pq <- apply(stage.year.EYT$snp, 2, function(x) var(x))     
  EYT.Genic <- sapply(EYT.alpha, function(x) EYT.pq %*% x^2)   # 2pqα^2
  EYT.Genic.var <- EYT.Genic*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)
  EYT.randGenic.true <- genParam(EYT.year)$randGenicVarG
  
  # EYT1 - 10 x 10
  snp.EYT1 <- stage.year.EYT$snp[1:10,]
  vector.year <- rep(16:21, each=10)
  EYT1.year <- EYT1[vector.year == year]
  EYT1.true.Genetic <- varG(EYT1.year)
  EYT1.true.Genic <- genicVarG(EYT1.year)
  EYT1.Genetic <- lapply(sample.effect, function(x) snp.EYT1 %*% eigenvectors %*% x)  # Wvα*
  EYT1.Genetic.est <- sapply(EYT1.Genetic, function(x) var(x))
  EYT1.Genetic.var <- EYT1.Genetic.est*(nrow(snp.EYT1)-1)/nrow(snp.EYT1)
  EYT1.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  EYT1.pq <- apply(snp.EYT1, 2, function(x) var(x))     
  EYT1.Genic <- sapply(EYT1.alpha, function(x) EYT1.pq %*% x^2)   # 2pqα^2
  EYT1.Genic.var <- EYT1.Genic*(nrow(snp.EYT1)-1)/nrow(snp.EYT1)
  EYT1.randGenic.true <- genParam(EYT1.year)$randGenicVarG
  
  # EYT2 - 10 x 10
  snp.EYT2 <- stage.year.EYT$snp[11:20,]
  vector.year <- rep(16:21, each=10)
  EYT2.year <- EYT2[vector.year == year]
  EYT2.true.Genetic <- varG(EYT2.year)
  EYT2.true.Genic <- genicVarG(EYT2.year)
  EYT2.Genetic <- lapply(sample.effect, function(x) snp.EYT2 %*% eigenvectors %*% x)  # Wvα*
  EYT2.Genetic.est <- sapply(EYT2.Genetic, function(x) var(x))
  EYT2.Genetic.var <- EYT2.Genetic.est*(nrow(snp.EYT2)-1)/nrow(snp.EYT2)
  EYT2.alpha <- lapply(sample.effect, function(x) eigenvectors %*% x)   # vα*
  EYT2.pq <- apply(snp.EYT2, 2, function(x) var(x))     
  EYT2.Genic <- sapply(EYT2.alpha, function(x) EYT2.pq %*% x^2)   # 2pqα^2
  EYT2.Genic.var <- EYT2.Genic*(nrow(snp.EYT2)-1)/nrow(snp.EYT2)
  EYT2.randGenic.true <- genParam(EYT2.year)$randGenicVarG
  
  save(stage.year.Parents, Parents.Genetic.var, Parents.Genic.var, Parents.randGenic.true,
       stage.year.F1, F1.Genetic.var, F1.Genic.var, F1.randGenic.true,
       stage.year.HDRW, HDRW.Genetic.var, HDRW.Genic.var, HDRW.randGenic.true,
       stage.year.PYT, PYT.Genetic.var, PYT.Genic.var, PYT.randGenic.true,
       stage.year.AYT, AYT.Genetic.var, AYT.Genic.var, AYT.randGenic.true,
       stage.year.EYT, EYT.Genetic.var, EYT.Genic.var, EYT.randGenic.true,
       EYT1.true.Genetic, EYT1.true.Genic, EYT1.Genetic.var, EYT1.Genic.var, EYT1.randGenic.true,
       EYT2.true.Genetic, EYT2.true.Genic, EYT2.Genetic.var, EYT2.Genic.var, EYT2.randGenic.true,
       file = paste0("./1000samples.10svd_year",year,".RData"))
}


##############################################################################################################
##############################################################################################################

