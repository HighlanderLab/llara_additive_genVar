
library(AlphaSimR)


source("./input_sim/GlobalParameters.R")
source("./input_sim/CreateParents.R")
source("./input_sim/FillPipeline.R")
output = data.frame(year = 1:40,
                    mean = numeric(40),
                    var = numeric(40),
                    genicVar = numeric(40),
                    accDH = numeric(40))
output.eff = data.frame(year = rep(16:21, each = 6), 
                        stage = rep(c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"),6),
                        mean = numeric(36),
                        var = numeric(36),
                        genicVar = numeric(36),
                        accDH = numeric(36))
parents.c <- F1c <- HDRWc <- PYTc <- AYTc <- EYTc <- NULL
i=1

# Cycle years 
for(year in 1:20){ 
  source("./input_sim/AdvanceYear.R") 
  
  # Crossing
  Parents = c(EYT2,EYT1,AYT)
  F1 = randCross(Parents,nCrosses)
  
  # Report results
  output$mean[year] = meanG(PYT)
  output$var[year] = varG(PYT)
  output$genicVar[year] = genicVarG(PYT)
  
  # Create training population
  if(year==16){
    trainPop = c(PYT,AYT,EYT1,EYT2)
    parents.c <- c(Parents)
    F1c <- c(F1)
    HDRWc <- c(HDRW)
    PYTc <- c(PYT)
    AYTc <- c(AYT)
    EYTc <- c(EYT1, EYT2)
  }else if(year>16){
    trainPop = c(trainPop,PYT,AYT,EYT1,EYT2)
    parents.c <- c(parents.c, Parents)
    F1c <- c(F1c, F1)
    HDRWc <- c(HDRWc, HDRW)
    PYTc <- c(PYTc, PYT)
    AYTc <- c(AYTc, AYT)
    EYTc <- c(EYTc, EYT1, EYT2)
  }
  
  # Report results per year and stages
  if(year > 15){
    EYT.year <- c(EYT1, EYT2)
    output.eff$mean[i:(i+5)] = c(meanG(Parents), meanG(F1), meanG(HDRW),
                                 meanG(PYT), meanG(AYT), meanG(EYT.year))
    output.eff$var[i:(i+5)] = c(varG(Parents), varG(F1), varG(HDRW),
                                varG(PYT), varG(AYT), varG(EYT.year))
    output.eff$genicVar[i:(i+5)] = c(genicVarG(Parents), genicVarG(F1), genicVarG(HDRW),
                                     genicVarG(PYT), genicVarG(AYT), genicVarG(EYT.year))
    i=i+6
    
    stage.year.Parents <- list(mean = meanG(Parents),
                               var = varG(Parents),
                               genicG = genicVarG(Parents),
                               id = Parents@id,
                               pheno = Parents@pheno,
                               fixed = Parents@fixEff,
                               snp = pullSnpGeno(Parents))
    stage.year.F1 <- list(mean = meanG(F1),
                          var = varG(F1),
                          genicG = genicVarG(F1),
                          id = F1@id,
                          pheno = F1@pheno,
                          fixed = F1@fixEff,
                          snp = pullSnpGeno(F1))
    stage.year.HDRW <- list(mean = meanG(HDRW),
                            var = varG(HDRW),
                            genicG = genicVarG(HDRW),
                            id = HDRW@id,
                            pheno = HDRW@pheno,
                            fixed = HDRW@fixEff,
                            snp = pullSnpGeno(HDRW))
    stage.year.PYT <- list(mean = meanG(PYT),
                           var = varG(PYT),
                           genicG = genicVarG(PYT),
                           id = PYT@id,
                           pheno = PYT@pheno,
                           fixed = PYT@fixEff,
                           snp = pullSnpGeno(PYT))
    stage.year.AYT <- list(mean = meanG(AYT),
                           var = varG(AYT),
                           genicG = genicVarG(AYT),
                           id = AYT@id,
                           pheno = AYT@pheno,
                           fixed = AYT@fixEff,
                           snp = pullSnpGeno(AYT))
    stage.year.EYT <- list(mean = meanG(EYT.year),
                           var = varG(EYT.year),
                           genicG = genicVarG(EYT.year),
                           id = EYT.year@id,
                           pheno = EYT.year@pheno,
                           fixed = EYT.year@fixEff,
                           snp = pullSnpGeno(EYT.year))
    save(stage.year.Parents, stage.year.F1, stage.year.HDRW, stage.year.PYT, stage.year.AYT,
         stage.year.EYT, file = paste0("./output_sim/stages_year", year, ".RData"))
  }
}

save.image("Burnin.RData")

