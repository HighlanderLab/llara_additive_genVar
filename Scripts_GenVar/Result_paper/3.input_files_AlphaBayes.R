
########################################################
#            Files input for AlphaBayes                #
#         Letícia Lara - January 30, 2020              #
########################################################


### Loading and manipulating the data
load("all_stages_and_year.RData")
trainSIM$fixed <- as.character(trainSIM$fixed)
trainSIM$id <- as.numeric(as.character(trainSIM$id))
head(trainSIM); str(trainSIM)
dim(snp);snp[1:5,1:5]


### PhenotypeTrainFile
pheno <- trainSIM[,c(1,2)]; head(pheno)
write.table(pheno, file = "PhenotypeTrain.txt", col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


### GenotypeTrainFile
rownames(snp) <- trainSIM$id
write.table(snp, file = "GenotypeTrain.txt", col.names = FALSE, quote = FALSE, sep = "\t")


### FixedEffectTrainFile
X <- model.matrix(~ 1 + fixed, data = trainSIM)
fixed <- cbind(trainSIM$id, data.frame(X[,-1])); head(fixed)
write.table(fixed, file = "HarvestTrain.txt", col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


