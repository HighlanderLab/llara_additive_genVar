

#######################################################################
#          Genetic and genic variances obtained from julia            #
#             Letícia Lara, Ivan Pocrnic, Gregor Gorjanc              #
#                             June, 2021                              #
#######################################################################

rm(list = ls())
library(AlphaSimR)             # Version 1.0.1
library(tidyverse)             # Version 1.3.0
library(dplyr)                 # Version 1.0.2


#### Loading and manipulating the data ####
QTLs = c(10, 100, 1000)
nSNPs = 500
scenGxY = c("noGxY", "GxY")
dir = "alphabayes"
startYear = 16
Yrs = 6

#### Running all repetitions ####
for(r in 1:10){
  for(Scen in scenGxY){
    for(nQTL in QTLs){
      load(paste0("./AlphaSimR/sim", Sim, "/Burnin_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
      load(paste0("./AlphaSimR/sim", Sim, "/all_stages_and_years_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
      
      full.sample <- matrix(as.numeric(scan(paste0("./", dir, "/", Scen, "/qtl", nQTL, "/sim", Sim, "/Samples_", dir, "_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".txt"),
                                            what = "numeric")), nrow = 10500)
      full.sample <- t(full.sample)
      dim(full.sample)
      sample.effect <- lapply(seq_len(nrow(full.sample)), function(i) full.sample[i,])
      length(sample.effect)
      
      
      #### Predictions of each stage and year of the breeding program ####
      vector.EYT <- rep(c(rep(1, 10), rep(2, 10)), Yrs)
      EYT1 <- EYTc[vector.EYT == 1]
      EYT2 <- EYTc[vector.EYT == 2]
      
      
      for(year in startYear:21){
        load(paste0("./AlphaSimR/sim", Sim, "/stages_year", year, "_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
        
        # Parents - 70 x 10500
        vector.year <- rep(startYear:21, each = 70)
        Parents.year <- parents.c[vector.year == year]
        Parents.true.Genetic <- varG(Parents.year)
        Parents.true.Genic <- genicVarG(Parents.year)
        Parents.Genetic <- lapply(sample.effect, function(x) stage.year.Parents$snp %*% x)    # Wα
        Parents.Genetic.est <- sapply(Parents.Genetic, function(x) var(x))
        Parents.Genetic.var <- Parents.Genetic.est*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)
        Parents.pq <- apply(stage.year.Parents$snp, 2, function(x) var(x))     
        Parents.Genic <- sapply(sample.effect, function(x) Parents.pq %*% x^2)   # var(g)α^2
        Parents.Genic.var <- Parents.Genic*(nrow(stage.year.Parents$snp)-1)/nrow(stage.year.Parents$snp)
        Parents.CovG_HW.true = genParam(Parents.year)$covG_HW  
        
        # F1 - 100 x 10500
        vector.year <- rep(startYear:21, each = 100)
        F1.year <- F1c[vector.year == year]
        F1.true.Genetic <- varG(F1.year)
        F1.true.Genic <- genicVarG(F1.year)
        F1.Genetic <- lapply(sample.effect, function(x) stage.year.F1$snp %*% x)  # Wα 
        F1.Genetic.est <- sapply(F1.Genetic, function(x) var(x))
        F1.Genetic.var <- F1.Genetic.est*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)
        F1.pq <- apply(stage.year.F1$snp, 2, function(x) var(x))     
        F1.Genic <-  sapply(sample.effect, function(x) F1.pq %*% x^2)   # var(g)α^2
        F1.Genic.var <- F1.Genic*(nrow(stage.year.F1$snp)-1)/nrow(stage.year.F1$snp)
        F1.CovG_HW.true = genParam(F1.year)$covG_HW  
        
        # HDRW - 10000 x 10500 from DH and 500 x 10500 from GS
        vector.year <- c(rep(startYear:20, each = 10000), rep(21, 500))
        HDRW.year <- HDRWc[vector.year == year]
        HDRW.true.Genetic <- varG(HDRW.year)
        HDRW.true.Genic <- genicVarG(HDRW.year)
        HDRW.Genetic <- lapply(sample.effect, function(x) stage.year.HDRW$snp %*% x)  # Wα
        HDRW.Genetic.est <- sapply(HDRW.Genetic, function(x) var(x))
        HDRW.Genetic.var <- HDRW.Genetic.est*(nrow(stage.year.HDRW$snp)-1)/nrow(stage.year.HDRW$snp)
        HDRW.pq <- apply(stage.year.HDRW$snp, 2, function(x) var(x))     
        HDRW.Genic <- sapply(sample.effect, function(x) HDRW.pq %*% x^2)   # var(g)α^2
        HDRW.Genic.var <- HDRW.Genic*(nrow(stage.year.HDRW$snp)-1)/nrow(stage.year.HDRW$snp)
        HDRW.CovG_HW.true = genParam(HDRW.year)$covG_HW  
        
        # PYT - 500 x 10500
        vector.year <- rep(startYear:21, each = 500)
        PYT.year <- PYTc[vector.year == year]
        PYT.true.Genetic <- varG(PYT.year)
        PYT.true.Genic <- genicVarG(PYT.year)
        PYT.Genetic <- lapply(sample.effect, function(x) stage.year.PYT$snp %*% x)  # Wα
        PYT.Genetic.est <- sapply(PYT.Genetic, function(x) var(x))
        PYT.Genetic.var <- PYT.Genetic.est*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)
        PYT.pq <- apply(stage.year.PYT$snp, 2, function(x) var(x))     
        PYT.Genic <- sapply(sample.effect, function(x) PYT.pq %*% x^2)   # var(g)α^2
        PYT.Genic.var <- PYT.Genic*(nrow(stage.year.PYT$snp)-1)/nrow(stage.year.PYT$snp)
        PYT.CovG_HW.true = genParam(PYT.year)$covG_HW  
        
        # AYT - 50 x 10500
        vector.year <- rep(startYear:21, each = 50)
        AYT.year <- AYTc[vector.year == year]
        AYT.true.Genetic <- varG(AYT.year)
        AYT.true.Genic <- genicVarG(AYT.year)
        AYT.Genetic <- lapply(sample.effect, function(x) stage.year.AYT$snp %*% x)  # Wα
        AYT.Genetic.est <- sapply(AYT.Genetic, function(x) var(x))
        AYT.Genetic.var <- AYT.Genetic.est*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)
        AYT.pq <- apply(stage.year.AYT$snp, 2, function(x) var(x))     
        AYT.Genic <- sapply(sample.effect, function(x) AYT.pq %*% x^2)   # var(g)α^2
        AYT.Genic.var <- AYT.Genic*(nrow(stage.year.AYT$snp)-1)/nrow(stage.year.AYT$snp)
        AYT.CovG_HW.true = genParam(AYT.year)$covG_HW  
        
        # EYT - 20 x 10500
        vector.year <- rep(startYear:21, each = 20)
        EYT.year <- EYTc[vector.year == year]
        EYT.true.Genetic <- varG(EYT.year)
        EYT.true.Genic <- genicVarG(EYT.year)
        EYT.Genetic <- lapply(sample.effect, function(x) stage.year.EYT$snp %*% x)  # Wα
        EYT.Genetic.est <- sapply(EYT.Genetic, function(x) var(x))
        EYT.Genetic.var <- EYT.Genetic.est*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)
        EYT.pq <- apply(stage.year.EYT$snp, 2, function(x) var(x))     
        EYT.Genic <- sapply(sample.effect, function(x) EYT.pq %*% x^2)   # var(g)α^2
        EYT.Genic.var <- EYT.Genic*(nrow(stage.year.EYT$snp)-1)/nrow(stage.year.EYT$snp)
        EYT.CovG_HW.true = genParam(EYT.year)$covG_HW  
        
        # EYT1 - 10 x 10500
        snp.EYT1 <- stage.year.EYT$snp[1:10,]
        vector.year <- rep(startYear:21, each = 10)
        EYT1.year <- EYT1[vector.year == year]
        gp_EYT1.year = genParam(EYT1.year)
        stage.year.EYT1 = list(mean = meanG(EYT1.year),
                               var = varG(EYT1.year),
                               genicG = genicVarG(EYT1.year),
                               genicA = genicVarA(EYT1.year),
                               covG_HW = gp_EYT1.year$covG_HW,
                               id = EYT1.year@id,
                               pheno = EYT1.year@pheno,
                               fixed = EYT1.year@fixEff,
                               snp = pullSnpGeno(EYT1.year))
        EYT1.true.Genetic <- varG(EYT1.year)
        EYT1.true.Genic <- genicVarG(EYT1.year)
        EYT1.Genetic <- lapply(sample.effect, function(x) stage.year.EYT1$snp %*% x)  # Wα
        EYT1.Genetic.est <- sapply(EYT1.Genetic, function(x) var(x))
        EYT1.Genetic.var <- EYT1.Genetic.est*(nrow(stage.year.EYT1$snp)-1)/nrow(stage.year.EYT1$snp)
        EYT1.pq <- apply(stage.year.EYT1$snp, 2, function(x) var(x))     
        EYT1.Genic <- sapply(sample.effect, function(x) EYT1.pq %*% x^2)   # var(g)α^2
        EYT1.Genic.var <- EYT1.Genic*(nrow(stage.year.EYT1$snp)-1)/nrow(stage.year.EYT1$snp)
        EYT1.CovG_HW.true = genParam(EYT1.year)$covG_HW  
        
        # EYT2 - 10 x 10500
        snp.EYT2 <- stage.year.EYT$snp[11:20,]
        vector.year <- rep(startYear:21, each = 10)
        EYT2.year <- EYT2[vector.year == year]
        gp_EYT2.year = genParam(EYT2.year)
        stage.year.EYT2 = list(mean = meanG(EYT2.year),
                               var = varG(EYT2.year),
                               genicG = genicVarG(EYT2.year),
                               genicA = genicVarA(EYT2.year),
                               covG_HW = gp_EYT2.year$covG_HW,
                               id = EYT2.year@id,
                               pheno = EYT2.year@pheno,
                               fixed = EYT2.year@fixEff,
                               snp = pullSnpGeno(EYT2.year))  
        EYT2.true.Genetic <- varG(EYT2.year)
        EYT2.true.Genic <- genicVarG(EYT2.year)
        EYT2.Genetic <- lapply(sample.effect, function(x) stage.year.EYT2$snp %*% x)  # Wα
        EYT2.Genetic.est <- sapply(EYT2.Genetic, function(x) var(x))
        EYT2.Genetic.var <- EYT2.Genetic.est*(nrow(stage.year.EYT2$snp)-1)/nrow(stage.year.EYT2$snp)
        EYT2.pq <- apply(stage.year.EYT2$snp, 2, function(x) var(x))     
        EYT2.Genic <- sapply(sample.effect, function(x) EYT2.pq %*% x^2)   # var(g)α^2
        EYT2.Genic.var <- EYT2.Genic*(nrow(stage.year.EYT2$snp)-1)/nrow(stage.year.EYT2$snp)
        EYT2.CovG_HW.true = genParam(EYT2.year)$covG_HW  
        
        save(Parents.true.Genetic, Parents.true.Genic, Parents.Genetic.var, Parents.Genic.var, Parents.CovG_HW.true, 
             F1.true.Genetic, F1.true.Genic, F1.Genetic.var, F1.Genic.var, F1.CovG_HW.true,
             HDRW.true.Genetic, HDRW.true.Genic, HDRW.Genetic.var, HDRW.Genic.var, HDRW.CovG_HW.true,
             PYT.true.Genetic, PYT.true.Genic, PYT.Genetic.var, PYT.Genic.var, PYT.CovG_HW.true,
             AYT.true.Genetic, AYT.true.Genic, AYT.Genetic.var, AYT.Genic.var, AYT.CovG_HW.true,
             EYT.true.Genetic, EYT.true.Genic, EYT.Genetic.var, EYT.Genic.var, EYT.CovG_HW.true, 
             EYT1.true.Genetic, EYT1.true.Genic, EYT1.Genetic.var, EYT1.Genic.var, EYT1.CovG_HW.true,
             EYT2.true.Genetic, EYT2.true.Genic, EYT2.Genetic.var, EYT2.Genic.var, EYT2.CovG_HW.true,
             file = paste0("./", dir, "/estimates/sim", Sim, "/Estimates_", dir, "_year", year, "_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
      }
    }
  }
}


