
sizeTrain = 16
Yrs = 6

if(intGxY == "noGxY"){
  pValueGxY = NULL
  varGE = 0
} else if(intGxY == "GxY"){
  pValueGxY = runif(5)        # generates random deviates from uniform distribution for FillPipeline.R script
  varGE = 0.2
}


print(paste0("GlobalParameters for REP ", r, ", int ", intGxY, ", and QTL ", nQTLs))
source("GlobalParameters.R")

print(paste0("CreateParents for REP ", r, ", int ", intGxY, ", and QTL ", nQTLs))
source("CreateParents.R")

print(paste0("FillPipeline for REP ", r, ", int ", intGxY, ", and QTL ", nQTLs))
source("FillPipeline.R")
output = data.frame(year = 1:40,
                    rep = rep(r, 40),
                    mean = numeric(40),
                    var = numeric(40),
                    genicVar = numeric(40),
                    covG_HW = numeric(40),
                    accDH = numeric(40))
output.eff = data.frame(year = rep(sizeTrain:21, each = 6), 
                        stage = rep(c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"), Yrs),
                        mean = numeric(Yrs*6),
                        var = numeric(Yrs*6),
                        genicVar = numeric(Yrs*6),
                        covG_HW = numeric(Yrs*6),
                        accDH = numeric(Yrs*6))
parents.c = F1c = HDRWc = PYTc = AYTc = EYTc = NULL
i = 1

# Cycle years 
if(intGxY == "noGxY"){
  pValueGxY = NULL
}else if(intGxY == "GxY"){
  pValueGxY = runif(21)        # generates random deviates from uniform distribution
}
print(paste0("AdvanceYear for REP", r, ", int ", intGxY, ", and QTL ", nQTLs))

for(year in 1:20){ 
  source("AdvanceYear.R") 
  
  # Crossing
  Parents = c(EYT2, EYT1, AYT)
  F1 = randCross(Parents, nCrosses)
  
  # Report results
  output$mean[year] = meanG(PYT)
  output$var[year] = varG(PYT)
  output$genicVar[year] = genicVarG(PYT)
  
  # Create training population
  if(year == sizeTrain){
    trainPop = c(PYT, AYT, EYT1, EYT2)
    parents.c = c(Parents)
    F1c = c(F1)
    HDRWc = c(HDRW)
    PYTc = c(PYT)
    AYTc = c(AYT)
    EYTc = c(EYT1, EYT2)
  }else if(year > sizeTrain){
    trainPop = c(trainPop, PYT, AYT, EYT1, EYT2)
    parents.c = c(parents.c, Parents)
    F1c = c(F1c, F1)
    HDRWc = c(HDRWc, HDRW)
    PYTc = c(PYTc, PYT)
    AYTc = c(AYTc, AYT)
    EYTc = c(EYTc, EYT1, EYT2)
  }
  
  # Report results per year and stages
  if(year >= sizeTrain){
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
  }
}

save.image(paste0("./AlphaSimR/sim", r, "/Burnin_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".RData"))


