
stage.year = list(Parents = list(), F1 = list(), HDRW = list(),
                  PYT = list(), AYT = list(), EYT = list()) 

# Cycle years to make more advanced parents
for(year in 21){               # Change to any number of desired years
  gsModel = RRBLUP(trainPop)
  print(paste0("AdvanceYearGS for REP ", r, ", int ", intGxY, ", and QTL ", nQTLs))
  source("AdvanceYearGS.R")    # Advances yield trials by a year
  
  # Crossing
  Parents = c(EYT2, EYT1, AYT)
  F1 = randCross(Parents, 100)

  # Report results
  output$mean[year] = meanG(PYT)
  output$var[year] = varG(PYT)
  output$genicVar[year] = genicVarG(PYT)

  EYT.year = c(EYT1, EYT2)
  gp_Parents = genParam(Parents)
  gp_F1 = genParam(F1)
  gp_HDRW = genParam(HDRW)
  gp_PYT = genParam(PYT)
  gp_AYT = genParam(AYT)
  gp_EYT.year = genParam(EYT.year)
  
  output.eff$mean[i:(i+5)] = c(meanG(Parents), meanG(F1), meanG(HDRW),
                               meanG(PYT), meanG(AYT), meanG(EYT.year))
  output.eff$var[i:(i+5)] = c(varG(Parents), varG(F1), varG(HDRW),
                              varG(PYT), varG(AYT), varG(EYT.year))
  output.eff$genicVar[i:(i+5)] = c(genicVarG(Parents), genicVarG(F1), genicVarG(HDRW),
                                   genicVarG(PYT), genicVarG(AYT), genicVarG(EYT.year))
  output.eff$covG_HW[i:(i+5)] = c(gp_Parents$covG_HW, gp_F1$covG_HW, gp_HDRW$covG_HW,
                                  gp_PYT$covG_HW, gp_AYT$covG_HW, gp_EYT.year$covG_HW)
  i = i + 6
  
  # Update training population
  trainPop = c(trainPop,PYT,AYT,EYT1,EYT2)
  parents.c = c(parents.c, Parents)
  F1c = c(F1c, F1)
  HDRWc = c(HDRWc, HDRW)
  PYTc = c(PYTc, PYT)
  AYTc = c(AYTc, AYT)
  EYTc = c(EYTc, EYT1, EYT2)

  # Saving files for each scenario and year
  stage.year.Parents = list(mean = meanG(Parents),
                            var = varG(Parents),
                            genicG = genicVarG(Parents),
                            genicA = genicVarA(Parents),
                            covG_HW = gp_Parents$covG_HW,
                            id = Parents@id,
                            pheno = Parents@pheno,
                            fixed = Parents@fixEff,
                            snp = pullSnpGeno(Parents))
  stage.year.F1 = list(mean = meanG(F1),
                       var = varG(F1),
                       genicG = genicVarG(F1),
                       genicA = genicVarA(F1),
                       covG_HW = gp_F1$covG_HW,
                       id = F1@id,
                       pheno = F1@pheno,
                       fixed = F1@fixEff,
                       snp = pullSnpGeno(F1))
  stage.year.HDRW = list(mean = meanG(HDRW),
                         var = varG(HDRW),
                         genicG = genicVarG(HDRW),
                         genicA = genicVarA(HDRW),
                         covG_HW = gp_HDRW$covG_HW,
                         id = HDRW@id,
                         pheno = HDRW@pheno,
                         fixed = HDRW@fixEff,
                         snp = pullSnpGeno(HDRW))
  stage.year.PYT = list(mean = meanG(PYT),
                        var = varG(PYT),
                        genicG = genicVarG(PYT),
                        genicA = genicVarA(PYT),
                        covG_HW = gp_PYT$covG_HW,
                        id = PYT@id,
                        pheno = PYT@pheno,
                        fixed = PYT@fixEff,
                        snp = pullSnpGeno(PYT))
  stage.year.AYT = list(mean = meanG(AYT),
                        var = varG(AYT),
                        genicG = genicVarG(AYT),
                        genicA = genicVarA(AYT),
                        covG_HW = gp_AYT$covG_HW,
                        id = AYT@id,
                        pheno = AYT@pheno,
                        fixed = AYT@fixEff,
                        snp = pullSnpGeno(AYT))
  stage.year.EYT = list(mean = meanG(EYT.year),
                        var = varG(EYT.year),
                        genicG = genicVarG(EYT.year),
                        genicA = genicVarA(EYT.year),
                        covG_HW = gp_EYT.year$covG_HW,
                        id = EYT.year@id,
                        pheno = EYT.year@pheno,
                        fixed = EYT.year@fixEff,
                        snp = pullSnpGeno(EYT.year))
  save(stage.year.Parents, stage.year.F1, stage.year.HDRW, stage.year.PYT, stage.year.AYT, stage.year.EYT, 
       file = paste0("./AlphaSimR/sim", r, "/stages_year", year, "_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".RData"))
  
  # Saving trainPop information for each stage of complete data (all years)
  gp_parents.c = genParam(parents.c)
  gp_F1c = genParam(F1c)
  gp_HDRWc = genParam(HDRWc)
  gp_PYTc = genParam(PYTc)
  gp_AYTc = genParam(AYTc)
  gp_EYTc = genParam(EYTc)
  
  stage.year$Parents = list(mean = meanG(parents.c),
                            var = varG(parents.c),
                            genicG = genicVarG(parents.c),
                            genicA = genicVarA(parents.c),
                            covG_HW = gp_parents.c$covG_HW,
                            id = parents.c@id,
                            pheno = parents.c@pheno,
                            fixed = parents.c@fixEff,
                            snp = pullSnpGeno(parents.c))
  stage.year$F1 = list(mean = meanG(F1c),
                       var = varG(F1c),
                       genicG = genicVarG(F1c),
                       genicA = genicVarA(F1c),
                       covG_HW = gp_F1c$covG_HW,
                       id = F1c@id,
                       pheno = F1c@pheno,
                       fixed = F1c@fixEff,
                       snp = pullSnpGeno(F1c))
  stage.year$HDRW = list(mean = meanG(HDRWc),
                         var = varG(HDRWc),
                         genicG = genicVarG(HDRWc),
                         genicA = genicVarA(HDRWc),
                         covG_HW = gp_HDRWc$covG_HW,
                         id = HDRWc@id,
                         pheno = HDRWc@pheno,
                         fixed = HDRWc@fixEff,
                         snp = pullSnpGeno(HDRWc))
  stage.year$PYT = list(mean = meanG(PYTc),
                        var = varG(PYTc),
                        genicG = genicVarG(PYTc),
                        genicA = genicVarA(PYTc),
                        covG_HW = gp_PYTc$covG_HW,
                        id = PYTc@id,
                        pheno = PYTc@pheno,
                        fixed = PYTc@fixEff,
                        snp = pullSnpGeno(PYTc))
  stage.year$AYT = list(mean = meanG(AYTc),
                        var = varG(AYTc),
                        genicG = genicVarG(AYTc),
                        genicA = genicVarA(AYTc),
                        covG_HW = gp_AYTc$covG_HW,
                        id = AYTc@id,
                        pheno = AYTc@pheno,
                        fixed = AYTc@fixEff,
                        snp = pullSnpGeno(AYTc))
  stage.year$EYT = list(mean = meanG(EYTc),
                        var = varG(EYTc),
                        genicG = genicVarG(EYTc),
                        genicA = genicVarA(EYTc),
                        covG_HW = gp_EYTc$covG_HW,
                        id = EYTc@id,
                        pheno = EYTc@pheno,
                        fixed = EYTc@fixEff,
                        snp = pullSnpGeno(EYTc))
}


# Getting information from train population 
trainSIM = data.frame("id" = trainPop@id,
                      "pheno" = trainPop@pheno,
                      "year"= trainPop@fixEff,
                      "stage"= rep(rep(c("PYT", "AYT", "EYT1", "EYT2"), c(nPYT, nAYT, nEYT, nEYT)), Yrs))
snp = pullSnpGeno(trainPop)

saveRDS(output, paste0("./AlphaSimR/sim", r, "/Base_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".rds"))

save(trainPop, trainSIM, snp, parents.c, F1c, HDRWc, PYTc, AYTc, EYTc, stage.year, output.eff,
     file = paste0("./AlphaSimR/sim", r, "/all_stages_and_years_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".RData"))

