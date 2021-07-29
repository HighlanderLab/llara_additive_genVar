
########################################################
#      Generating files for AlphaBayes using SVD       #
#     Letícia Lara, Ivan Pocrnic, Gregor Gorjanc       #
#                     March, 2020                      #
########################################################

library(AlphaSimR)           # Version 0.11.0
library(irlba)


#### Loading and manipulating the data ####
load("../all_stages_and_year.RData")

trainSIM$fixed <- as.character(trainSIM$fixed)
trainSIM$id <- as.numeric(as.character(trainSIM$id))
trainSIM$id2 <- 1:nrow(snp)
head(trainSIM); str(trainSIM)
dim(snp); snp[1:5,1:5]


#### Singular value decomposition ####
# irlba_10pc$v == irlba_10pc$rotation
# Svd.10PC_snp == irlba_10pc$x == mat_10pc

irlba_10pc <- prcomp_irlba(snp, n=10)
mat_10pc <- irlba_10pc$x
rownames(mat_10pc) <- trainSIM$id
write.table(mat_10pc, file = "GenotypeTrain10PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_10pc, file = "prcomp.irlba_10pc.RData")

irlba_50pc <- prcomp_irlba(snp, n=50)
mat_50pc <- irlba_50pc$x
rownames(mat_50pc) <- trainSIM$id
write.table(mat_50pc, file = "GenotypeTrain50PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_50pc, file = "prcomp.irlba_50pc.RData")

irlba_100pc <- prcomp_irlba(snp, n=100)
mat_100pc <- irlba_100pc$x
rownames(mat_100pc) <- trainSIM$id
write.table(mat_100pc, file = "GenotypeTrain100PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_100pc, file = "prcomp.irlba_100pc.RData")

irlba_500pc <- prcomp_irlba(snp, n=500)
mat_500pc <- irlba_500pc$x
rownames(mat_500pc) <- trainSIM$id
write.table(mat_500pc, file = "GenotypeTrain500PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_500pc, file = "prcomp.irlba_500pc.RData")

irlba_1000pc <- prcomp_irlba(snp, n=1000)
mat_1000pc <- irlba_1000pc$x
rownames(mat_1000pc) <- trainSIM$id
write.table(mat_1000pc, file = "GenotypeTrain1000PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_1000pc, file = "prcomp.irlba_1000pc.RData")

irlba_2000pc <- prcomp_irlba(snp, n=2000)
mat_2000pc <- irlba_2000pc$x
rownames(mat_2000pc) <- trainSIM$id
write.table(mat_2000pc, file = "GenotypeTrain2000PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_2000pc, file = "prcomp.irlba_2000pc.RData")

irlba_3410pc <- prcomp_irlba(snp, n=3410)
mat_3410pc <- irlba_3410pc$x
rownames(mat_3410pc) <- trainSIM$id
write.table(mat_3410pc, file = "GenotypeTrain3410PC.txt", col.names = FALSE, quote = FALSE, sep = "\t")
save(irlba_3410pc, file = "prcomp.irlba_3410pc.RData")


#### Reading files ####
mat_10pc <- read.table("./10pc/GenotypeTrain10PC.txt")
dim(mat_10pc); class(mat_10pc)
summary(as.vector(as.matrix(mat_10pc[,2:11])))

mat_50pc <- read.table("./50pc/GenotypeTrain50PC.txt")
dim(mat_50pc); class(mat_50pc)
summary(as.vector(as.matrix(mat_50pc[,2:51])))

mat_100pc <- read.table("./100pc/GenotypeTrain100PC.txt")
dim(mat_100pc); class(mat_100pc)
summary(as.vector(as.matrix(mat_100pc[,2:101])))

mat_500pc <- read.table("./500pc/GenotypeTrain500PC.txt")
dim(mat_500pc); class(mat_500pc)
summary(as.vector(as.matrix(mat_500pc[,2:501])))

mat_1000pc <- read.table("./1000pc/GenotypeTrain1000PC.txt")
dim(mat_1000pc); class(mat_1000pc)
summary(as.vector(as.matrix(mat_1000pc[,2:1001])))

mat_2000pc <- read.table("./2000pc/GenotypeTrain2000PC.txt")
dim(mat_2000pc); class(mat_2000pc)
summary(as.vector(as.matrix(mat_2000pc[,2:2001])))

mat_3410pc <- read.table("./3410pc/GenotypeTrain3410PC.txt")
dim(mat_3410pc); class(mat_3410pc)
summary(as.vector(as.matrix(mat_3410pc[,2:3411])))


##############################################################################################################
##############################################################################################################
