
########################################################
#            Files input for AlphaBayes                #
#             Letícia Lara - June, 2020                #
########################################################

library(AlphaSimR)

for(simulation in 1:10){
  print(paste0("This is loop ", simulation))
  REP = paste0("_Sim", simulation)
  
  ### Loading and manipulating the data
  load(paste0("./AlphaSimR/Results/all_stages_and_years", REP, ".RData"))
  trainSIM$fixed <- as.character(trainSIM$fixed)
  trainSIM$id <- as.numeric(as.character(trainSIM$id))

  ### PhenotypeTrainFile
  pheno <- trainSIM[,c(1,2)]
  write.table(pheno, file = paste0("./AlphaBayes/Simulation", simulation, "/PhenotypeTrain", REP, ".txt"), 
              col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

  ### GenotypeTrainFile
  rownames(snp) <- trainSIM$id
  write.table(snp, file = paste0("./AlphaBayes/Simulation", simulation, "/GenotypeTrain", REP, ".txt"), 
              col.names = FALSE, quote = FALSE, sep = "\t")

  ### FixedEffectTrainFile
  X <- model.matrix(~ 1 + fixed, data = trainSIM)
  fixed <- cbind(trainSIM$id, data.frame(X[,-1]))
  write.table(fixed, file = paste0("./AlphaBayes/Simulation", simulation, "/HarvestTrain", REP, ".txt"), 
              col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  
  ### SpecificationFile 
  sink(file = paste0("./AlphaBayes/Simulation", simulation, "/fullBayes/AlphaBayesSpec.txt"))
  cat(paste0("PhenotypeTrainFile         , PhenotypeTrain", REP, ".txt\n"))
  cat(paste0("FixedEffectTrainFile       , HarvestTrain", REP, ".txt\n"))
  cat(paste0("NumberOfFixedEffectLevels  , ", (ncol(fixed)-1), "\n"))
  cat(paste0("GenotypeTrainFile          , GenotypeTrain", REP, ".txt\n"))
  cat(paste0("NumberOfMarkers            , ", ncol(snp), "\n"))
  cat("GenotypeScalingMethod      , 1\n")
  cat("EstimateVariances          , Yes\n")
  cat("EstimationMethod           , RidgeSample\n")
  cat("NumberOfIterations         , 100000\n")
  cat("NumberOfBurnInIterations   , 10000\n")
  cat("Stop\n")
  sink()
}


#### After running AlphaBayes as fullBayes approach -> RUN empiricalBayes ####
nFixed = 5
nMarkers = 10500

for(simulation in 1:10){
  print(paste0("This is loop ", simulation))
  REP = paste0("_Sim", simulation)
  
  genVar = read.table(paste0("./AlphaBayes-Estimation/Simulation", simulation, "/fullBayes/MarkerVarianceEstimate.txt")) 
  genVar = round(as.numeric(genVar*nMarkers), 7)
  resVar = read.table(paste0("./AlphaBayes-Estimation/Simulation", simulation, "/fullBayes/ResidualVarianceEstimate.txt")) 
  resVar = round(as.numeric(resVar), 7)
  
  ### SpecificationFile 
  sink(file = paste0("./AlphaBayes-Estimation/Simulation", simulation, "/empiricalBayes/AlphaBayesSpec.txt"))
  cat(paste0("PhenotypeTrainFile         , PhenotypeTrain", REP, ".txt\n"))
  cat(paste0("FixedEffectTrainFile       , HarvestTrain", REP, ".txt\n"))
  cat(paste0("NumberOfFixedEffectLevels  , ", nFixed, "\n"))
  cat(paste0("GenotypeTrainFile          , GenotypeTrain", REP, ".txt\n"))
  cat(paste0("NumberOfMarkers            , ", nMarkers, "\n"))
  cat("GenotypeScalingMethod      , 1\n")
  cat("EstimateVariances          , No\n")
  cat("EstimationMethod           , RidgeSample\n")
  cat(paste0("GeneticVariance            , ", genVar, "\n"))
  cat(paste0("ResidualVariance           , ", resVar, "\n"))
  cat("NumberOfIterations         , 100000\n")
  cat("NumberOfBurnInIterations   , 10000\n")
  cat("Stop\n")
  sink()
}

