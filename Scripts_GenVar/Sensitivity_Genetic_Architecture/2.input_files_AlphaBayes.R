
########################################################
#            Files input for AlphaBayes                #
#            Letícia Lara - March, 2021                #
########################################################

library(AlphaSimR)             # Version 1.0.1

QTLs = c(10, 100, 1000)
nSNPs = 500
scenGxY = c("noGxY", "GxY")

for(r in 1:10){
  for(intGxY in scenGxY){
    for(nQTLs in QTLs){
      ### Loading and manipulating the data
      load(paste0("./AlphaSimR/sim", r, "/all_stages_and_years_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".RData"))
      trainSIM$year <- as.character(trainSIM$year)
      trainSIM$id <- as.numeric(as.character(trainSIM$id))

      ### PhenotypeTrainFile
      pheno <- trainSIM[,c(1,2)]
      write.table(pheno, file = paste0("./alphabayes/", intGxY, "/qtl", nQTLs, "/sim", r, "/Phenotype_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt"), 
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

      ### GenotypeTrainFile
      rownames(snp) <- trainSIM$id
      write.table(snp, file = paste0("./alphabayes/", intGxY, "/qtl", nQTLs, "/sim", r, "/Genotype_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt"), 
                  sep = "\t", col.names = FALSE, quote = FALSE)

      ### FixedEffectTrainFile
      X <- model.matrix(~ 1 + year, data = trainSIM)
      fixed <- cbind(trainSIM$id, data.frame(X[,-1]))
      write.table(fixed, file = paste0("./alphabayes/", intGxY, "/qtl", nQTLs, "/sim", r, "/Harvest_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt"),
                  col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

      ### SpecificationFile
sink(file = paste0("./alphabayes/", intGxY, "/qtl", nQTLs, "/sim", r, "/AlphaBayesSpec.txt"))
cat(paste0("PhenotypeTrainFile         , Phenotype_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt\n"))
cat(paste0("FixedEffectTrainFile       , Harvest_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt\n"))
cat(paste0("NumberOfFixedEffectLevels  , ", (ncol(fixed)-1), "\n"))
cat(paste0("GenotypeTrainFile          , Genotype_AlphaBayes_", intGxY, "_nQTLs", nQTLs, "_Sim", r, ".txt\n"))
cat(paste0("NumberOfMarkers            , ", ncol(snp), "\n"))
cat("GenotypeScalingMethod	       , 1\n")
cat("EstimateVariances                 , Yes\n")
cat("EstimationMethod                  , RidgeSample\n")
cat("NumberOfIterations                , 100000\n")
cat("NumberOfBurnInIterations          , 10000\n")
cat("Stop\n")
sink()
    }
  }
}





