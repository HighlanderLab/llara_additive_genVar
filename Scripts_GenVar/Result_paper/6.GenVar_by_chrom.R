
#######################################################################
#          Estimating genetic, genic, and LD by chromosome            #
#                 Simulated wheat breeding program                    #
#                  Leticia Lara - February, 2020                      #
#######################################################################

library(AlphaSimR)           # Version 0.11.0
library(tidyverse)
library(dplyr)
INLA:::inla.dynload.workaround()


#############################################
# Estimates for Trait effects - True values
#### True values for EYT ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/EYT_year21.RData")

vector.EYT <- rep(c(rep(1,10), rep(2,10)),6)
EYT1 <- EYTc[vector.EYT == 1]
vector.year <- rep(16:21, each=10)
EYT1.year <- EYT1[vector.year == 21]

byQTL.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  byQTL.chr <- NULL
  for(i in 1:100){
    byQTL <- pullQtlGeno(EYT1.year, chr=c)[,i] * SP$traits[[1]]@addEff[i+j]
    byQTL.matrix <- cbind(byQTL.matrix, byQTL)
    byQTL.chr <- cbind(byQTL.chr, byQTL)
  }
  matrix.chr[c] <- list(byQTL.chr)
  j=i+j
}
QTL.matrix <- popVar(byQTL.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
QTL.matrix[1:5,1:5]

vector <- rbind(seq(1, 2001, 100), seq(100, 2100, 100))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(QTL.matrix[i:j, i:j])
  genic.value <- sum(diag(QTL.matrix[i:j,i:j]))
  LD.without <- (sum(QTL.matrix[,i:j]) - sum(QTL.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all
varG(EYT1.year)
sum(unlist(chr.genic)); genicVarG(EYT1.year)

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_trait_true_EYT.RData")


#### True values for HDRW ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load(paste0("../output_sim/HDRW_year18.RData"))

vector.year <- c(rep(16:20, each=10000), rep(21, each=500))
HDRW.year <- HDRWc[vector.year == 18]

byQTL.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  byQTL.chr <- NULL
  for(i in 1:100){
    byQTL <- pullQtlGeno(HDRW.year, chr=c)[,i] * SP$traits[[1]]@addEff[i+j]
    byQTL.matrix <- cbind(byQTL.matrix, byQTL)
    byQTL.chr <- cbind(byQTL.chr, byQTL)
  }
  matrix.chr[c] <- list(byQTL.chr)
  j=i+j
}
QTL.matrix <- popVar(byQTL.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
QTL.matrix[1:5,1:5]

vector <- rbind(seq(1, 2001, 100), seq(100, 2100, 100))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(QTL.matrix[i:j, i:j])
  genic.value <- sum(diag(QTL.matrix[i:j,i:j]))
  LD.without <- (sum(QTL.matrix[,i:j]) - sum(QTL.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all
varG(HDRW.year)
sum(unlist(chr.genic)); genicVarG(HDRW.year)

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_trait_true_HDRW.RData")


#############################################
# Estimates for genome - QTLs
#### True values for EYT ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load(paste0("../output_sim/EYT_year21.RData"))

vector.EYT <- rep(c(rep(1,10), rep(2,10)),6)
EYT1 <- EYTc[vector.EYT == 1]
vector.year <- rep(16:21, each=10)
EYT1.year <- EYT1[vector.year == 21]

byQTL.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  byQTL.chr <- NULL
  for(i in 1:100){
    byQTL <- pullQtlGeno(EYT1.year, chr=c)[,i]
    byQTL.matrix <- cbind(byQTL.matrix, byQTL)
    byQTL.chr <- cbind(byQTL.chr, byQTL)
  }
  matrix.chr[c] <- list(byQTL.chr)
  j=i+j
}
QTL.matrix <- popVar(byQTL.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
QTL.matrix[1:5,1:5]

vector <- rbind(seq(1, 2001, 100), seq(100, 2100, 100))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(QTL.matrix[i:j, i:j])
  genic.value <- sum(diag(QTL.matrix[i:j,i:j]))
  LD.without <- (sum(QTL.matrix[,i:j]) - sum(QTL.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_genome_QTLs_EYT.RData")


#### True values for HDRW ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/HDRW_year18.RData")

vector.year <- c(rep(16:20, each=10000), rep(21, each=500))
HDRW.year <- HDRWc[vector.year == 18]

byQTL.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  byQTL.chr <- NULL
  for(i in 1:100){
    byQTL <- pullQtlGeno(HDRW.year, chr=c)[,i]
    byQTL.matrix <- cbind(byQTL.matrix, byQTL)
    byQTL.chr <- cbind(byQTL.chr, byQTL)
  }
  matrix.chr[c] <- list(byQTL.chr)
  j=i+j
}
QTL.matrix <- popVar(byQTL.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
QTL.matrix[1:5,1:5]

vector <- rbind(seq(1, 2001, 100), seq(100, 2100, 100))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(QTL.matrix[i:j, i:j])
  genic.value <- sum(diag(QTL.matrix[i:j,i:j]))
  LD.without <- (sum(QTL.matrix[,i:j]) - sum(QTL.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_genome_QTLs_HDRW.RData")


#############################################
# Estimates for Trait effects  - Est values
#### Estimated values for EYT ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/EYT_year21.RData")
load("../fullBayes/sample.effect_full.RData")

vector.EYT <- rep(c(rep(1,10), rep(2,10)),6)
EYT1 <- EYTc[vector.EYT == 1]
vector.year <- rep(16:21, each=10)
EYT1.year <- EYT1[vector.year == 21]

sample.genetic <- sample.genic <- NULL
sample.LD <- sample.LD.without <- NULL
for(k in 1:1000){
  print(paste0("#### Loop: ", k))
  bySNP.matrix <- NULL; j=0
  matrix.chr <- list()
  for(c in 1:21){
    bySNP.chr <- NULL
    for(i in 1:500){
      bySNP <- pullSnpGeno(EYT1.year, chr=c)[,i] * sample.effect[[k]][i+j]
      bySNP.matrix <- cbind(bySNP.matrix, bySNP)
      bySNP.chr <- cbind(bySNP.chr, bySNP)
    }
    matrix.chr[c] <- list(bySNP.chr)
    j=i+j
  }
  SNP.matrix <- popVar(bySNP.matrix)
  matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
  dim(SNP.matrix); SNP.matrix[1:5,1:5]
  
  vector <- rbind(seq(1, 10001, 500), seq(500, 10500, 500))
  chr.genetic <- chr.genic <- chr.LD.without <- NULL
  for(c in 1:21){
    i = vector[1,c]
    j = vector[2,c]
    genetic.value <- sum(SNP.matrix[i:j, i:j])
    genic.value <- sum(diag(SNP.matrix[i:j,i:j]))
    LD.without <- (sum(SNP.matrix[,i:j]) - sum(SNP.matrix[i:j,i:j]))
    chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
    chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
    chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
  }
  chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
  
  sample.genetic <- cbind(sample.genetic, unlist(chr.genetic)) 
  sample.genic <- cbind(sample.genic, unlist(chr.genic))
  sample.LD <- cbind(sample.LD, unlist(chr.LD))
  sample.LD.without <- cbind(sample.LD.without, unlist(chr.LD.without))
}
EYT.genetic <- rowMeans(sample.genetic)
EYT.genic <- rowMeans(sample.genic)
EYT.LD <- rowMeans(sample.LD)
EYT.LD.without <- rowMeans(sample.LD.without)
sum(c(EYT.genic, EYT.LD, EYT.LD.without))         # Number 0 - Genetic for all

save(sample.genetic, sample.genic, sample.LD, sample.LD.without,
     EYT.genetic, EYT.genic, EYT.LD, EYT.LD.without, file = "byChrom_trait_est_EYT.RData")


#### Estimated values for HDRW ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/HDRW_year18.RData")
load("../fullBayes/sample.effect_full.RData")

vector.year <- c(rep(16:20, each=10000), rep(21, each=500))
HDRW.year <- HDRWc[vector.year == 18]

sample.genetic <- sample.genic <- NULL
sample.LD <- sample.LD.without <- NULL
for(k in 1:1000){
  print(paste0("#### Loop: ", k))
  bySNP.matrix <- NULL; j=0
  matrix.chr <- list()
  for(c in 1:21){
    bySNP.chr <- NULL
    for(i in 1:500){
      bySNP <- pullSnpGeno(HDRW.year, chr=c)[,i] * sample.effect[[k]][i+j]
      bySNP.matrix <- cbind(bySNP.matrix, bySNP)
      bySNP.chr <- cbind(bySNP.chr, bySNP)
    }
    matrix.chr[c] <- list(bySNP.chr)
    j=i+j
  }
  SNP.matrix <- popVar(bySNP.matrix)
  matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
  dim(SNP.matrix); SNP.matrix[1:5,1:5]
  
  vector <- rbind(seq(1, 10001, 500), seq(500, 10500, 500))
  chr.genetic <- chr.genic <- chr.LD.without <- NULL
  for(c in 1:21){
    i = vector[1,c]
    j = vector[2,c]
    genetic.value <- sum(SNP.matrix[i:j, i:j])
    genic.value <- sum(diag(SNP.matrix[i:j,i:j]))
    LD.without <- (sum(SNP.matrix[,i:j]) - sum(SNP.matrix[i:j,i:j]))
    chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
    chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
    chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
  }
  chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
  
  sample.genetic <- cbind(sample.genetic, unlist(chr.genetic)) 
  sample.genic <- cbind(sample.genic, unlist(chr.genic))
  sample.LD <- cbind(sample.LD, unlist(chr.LD))
  sample.LD.without <- cbind(sample.LD.without, unlist(chr.LD.without))
}
HDRW.genetic <- rowMeans(sample.genetic)
HDRW.genic <- rowMeans(sample.genic)
HDRW.LD <- rowMeans(sample.LD)
HDRW.LD.without <- rowMeans(sample.LD.without)
sum(c(HDRW.genic, HDRW.LD, HDRW.LD.without))      # Number 0 - Genetic for all

save(sample.genetic, sample.genic, sample.LD, sample.LD.without,
     HDRW.genetic, HDRW.genic, HDRW.LD, HDRW.LD.without, file = "byChrom_trait_est_HDRW.RData")


#############################################
# Estimates for genome - SNPs
#### Estimated values for EYT ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/EYT_year21.RData")
load("../fullBayes/sample.effect_full.RData")

vector.EYT <- rep(c(rep(1,10), rep(2,10)),6)
EYT1 <- EYTc[vector.EYT == 1]
vector.year <- rep(16:21, each=10)
EYT1.year <- EYT1[vector.year == 21]

bySNP.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  bySNP.chr <- NULL
  for(i in 1:500){
    bySNP <- pullSnpGeno(EYT1.year, chr=c)[,i]
    bySNP.matrix <- cbind(bySNP.matrix, bySNP)
    bySNP.chr <- cbind(bySNP.chr, bySNP)
  }
  matrix.chr[c] <- list(bySNP.chr)
  j=i+j
}
SNP.matrix <- popVar(bySNP.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
dim(SNP.matrix); SNP.matrix[1:5,1:5]
  
vector <- rbind(seq(1, 10001, 500), seq(500, 10500, 500))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(SNP.matrix[i:j, i:j])
  genic.value <- sum(diag(SNP.matrix[i:j,i:j]))
  LD.without <- (sum(SNP.matrix[,i:j]) - sum(SNP.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_genome_SNPs_EYT.RData")


#### Estimated values for HDRW ####
rm(list=ls())
load("../Burnin.RData")
load("../all_stages_and_year.RData")
load("../output_sim/HDRW_year18.RData")
load("../fullBayes/sample.effect_full.RData")

vector.year <- c(rep(16:20, each=10000), rep(21, each=500))
HDRW.year <- HDRWc[vector.year == 18]

bySNP.matrix <- NULL; j=0
matrix.chr <- list()
for(c in 1:21){
  bySNP.chr <- NULL
  for(i in 1:500){
    bySNP <- pullSnpGeno(HDRW.year, chr=c)[,i]
    bySNP.matrix <- cbind(bySNP.matrix, bySNP)
    bySNP.chr <- cbind(bySNP.chr, bySNP)
  }
  matrix.chr[c] <- list(bySNP.chr)
  j=i+j
}
SNP.matrix <- popVar(bySNP.matrix)
matrix.chr <- lapply(matrix.chr, function(x) popVar(x))
dim(SNP.matrix); SNP.matrix[1:5,1:5]
  
vector <- rbind(seq(1, 10001, 500), seq(500, 10500, 500))
chr.genetic <- chr.genic <- chr.LD.without <- NULL
for(c in 1:21){
  i = vector[1,c]
  j = vector[2,c]
  genetic.value <- sum(SNP.matrix[i:j, i:j])
  genic.value <- sum(diag(SNP.matrix[i:j,i:j]))
  LD.without <- (sum(SNP.matrix[,i:j]) - sum(SNP.matrix[i:j,i:j]))
  chr.genetic[c] <- list(genetic.value)         # Number 1 - genetic estimates by chromosome
  chr.genic[c] <- list(genic.value)             # Number 2 - genic estimates by chromosome / Sum of chr.genic is genic for all chromosomes (E2)
  chr.LD.without[c] <- list(LD.without)         # Number 4 -  LD without the chromosome
}

chr.LD <- lapply(matrix.chr, function(x) sum(x) - sum(diag(x)))     # Number 3 - LD by chromosome
sum(c(unlist(chr.genic), unlist(chr.LD), unlist(chr.LD.without)))   # Number 0 - Genetic for all

save(chr.genetic, chr.genic, chr.LD, chr.LD.without, file = "byChrom_genome_SNPs_HDRW.RData")



##############################################################################################################
##############################################################################################################

########################################################
#### Estimates for trait effects - Estimated values ####
# Loading results from Eddie
rm(list=ls())
load("sample1.100_HDRW.RData")  
load("sample101.200_HDRW.RData")
load("sample201.300_HDRW.RData")
load("sample301.400_HDRW.RData")
load("sample401.500_HDRW.RData")
load("sample501.600_HDRW.RData")
load("sample601.700_HDRW.RData")
load("sample701.800_HDRW.RData")
load("sample801.900_HDRW.RData")
load("sample901.1000_HDRW.RData")

dim(sample1.genetic); dim(sample1.genic)
dim(sample1.LD); dim(sample1.LD.without)

sample.genetic <- cbind(sample1.genetic, sample2.genetic, sample3.genetic, sample4.genetic, 
                        sample5.genetic, sample6.genetic, sample7.genetic, sample8.genetic, 
                        sample9.genetic, sample10.genetic)
sample.genic <- cbind(sample1.genic, sample2.genic, sample3.genic, sample4.genic,
                      sample5.genic, sample6.genic, sample7.genic, sample8.genic,
                      sample9.genic, sample10.genic)
sample.LD <- cbind(sample1.LD, sample2.LD, sample3.LD, sample4.LD, sample5.LD,
                   sample6.LD, sample7.LD, sample8.LD, sample9.LD, sample10.LD)
sample.LD.without <- cbind(sample1.LD.without, sample2.LD.without, sample3.LD.without, sample4.LD.without,
                           sample5.LD.without, sample6.LD.without, sample7.LD.without, sample8.LD.without,
                           sample9.LD.without, sample10.LD.without)

HDRW.genetic <- rowMeans(sample.genetic)
HDRW.genic <- rowMeans(sample.genic)
HDRW.LD <- rowMeans(sample.LD)
HDRW.LD.without <- rowMeans(sample.LD.without)
sum(c(HDRW.genic, HDRW.LD, HDRW.LD.without))      # Number 0 - Genetic for all
sum(HDRW.genic)

save(sample.genetic, sample.genic, sample.LD, sample.LD.without,
     HDRW.genetic, HDRW.genic, HDRW.LD, HDRW.LD.without, file = "byChrom_trait_est_HDRW_fullBayes.RData")


##############################################################################################################
##############################################################################################################

#############################################
#### Graph - Trait effects ####
library(ggplot2)
library(tidyverse)
rm(list=ls())

load("byChrom_trait_true_EYT.RData")
EYT.genetic_trait_true <- unlist(chr.genetic)
EYT.genic_trait_true <- unlist(chr.genic)
EYT.LD_trait_true <- unlist(chr.LD)
EYT.LD.without_trait_true <- unlist(chr.LD.without)

load("byChrom_trait_true_HDRW.RData")
HDRW.genetic_trait_true <- unlist(chr.genetic)
HDRW.genic_trait_true <- unlist(chr.genic)
HDRW.LD_trait_true <- unlist(chr.LD)
HDRW.LD.without_trait_true <- unlist(chr.LD.without)

load("byChrom_trait_est_EYT.RData")
EYT.genetic_trait_est <- EYT.genetic
EYT.genic_trait_est <- EYT.genic
EYT.LD_trait_est <- EYT.LD
EYT.LD.without_trait_est <- EYT.LD.without
sample.EYT.genetic <- sample.genetic
sample.EYT.genic <- sample.genic
sample.EYT.LD <- sample.LD
sample.EYT.LD.without <- sample.LD.without

load("byChrom_trait_est_HDRW_fullBayes.RData")
HDRW.genetic_trait_est <- HDRW.genetic
HDRW.genic_trait_est <- HDRW.genic
HDRW.LD_trait_est <- HDRW.LD
HDRW.LD.without_trait_est <- HDRW.LD.without
sample.HDRW.genetic <- sample.genetic
sample.HDRW.genic <- sample.genic
sample.HDRW.LD <- sample.LD
sample.HDRW.LD.without <- sample.LD.without

rm(chr.genetic, chr.genic, chr.LD, chr.LD.without,
   sample.genetic, sample.genic, sample.LD, sample.LD.without,
   EYT.genetic, EYT.genic, EYT.LD, EYT.LD.without,
   HDRW.genetic, HDRW.genic, HDRW.LD, HDRW.LD.without)

length(EYT.genetic_trait_true)
dim(sample.EYT.genetic)


# Using total values
sum_sample.EYT.genetic <- colSums(rbind(sample.EYT.genic, sample.EYT.LD, sample.EYT.LD.without))
sum_sample.EYT.genic <- colSums(sample.EYT.genic)
sum_sample.EYT.LD <- colSums(sample.EYT.LD)
sum_sample.EYT.LD.without <- colSums(sample.EYT.LD.without)
sum_sample.HDRW.genetic <- colSums(rbind(sample.HDRW.genic, sample.HDRW.LD, sample.HDRW.LD.without))
sum_sample.HDRW.genic <- colSums(sample.HDRW.genic)
sum_sample.HDRW.LD <- colSums(sample.HDRW.LD)
sum_sample.HDRW.LD.without <- colSums(sample.HDRW.LD.without)

sum_true.EYT.genetic <- sum(EYT.genic_trait_true, EYT.LD_trait_true, EYT.LD.without_trait_true)
sum_true.EYT.genic <- sum(EYT.genic_trait_true)
sum_true.EYT.LD <- sum(EYT.LD_trait_true)
sum_true.EYT.LD.without <- sum(EYT.LD.without_trait_true)
sum_true.HDRW.genetic <- sum(HDRW.genic_trait_true, HDRW.LD_trait_true, HDRW.LD.without_trait_true)
sum_true.HDRW.genic <- sum(HDRW.genic_trait_true)
sum_true.HDRW.LD <- sum(HDRW.LD_trait_true)
sum_true.HDRW.LD.without <- sum(HDRW.LD.without_trait_true)

sum_true <- rep(c(sum_true.HDRW.genetic, sum_true.EYT.genetic,
                  sum_true.HDRW.genic, sum_true.EYT.genic,
                  sum_true.HDRW.LD, sum_true.EYT.LD,
                  sum_true.HDRW.LD.without, sum_true.EYT.LD.without), each=1000)

optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_roslin <- c("#5FA83B", "#B2E39C", "#C11773", "#DE99BF", "#3798D9", "#9DC4E0", "#DB760FFF", "#FAB493")

trait_est.graph <- data.frame("Estimates" = rep(c(rep("Genetic", 1000), rep("Genic", 1000),
                                                  rep("Within-LD", 1000), rep("Between-LD", 1000)), 2),
                              "Values" <- c(sum_sample.HDRW.genetic, sum_sample.HDRW.genic,
                                            sum_sample.HDRW.LD, sum_sample.HDRW.LD.without,
                                            sum_sample.EYT.genetic, sum_sample.EYT.genic, 
                                            sum_sample.EYT.LD, sum_sample.EYT.LD.without),
                              "Stage" <- c(rep("HDRW", 4000), rep("EYT", 4000)))
colnames(trait_est.graph) <- c("Estimates", "Values", "Stage")
trait_est.graph$Estimates <- factor(trait_est.graph$Estimates, 
                                    levels = c("Genetic", "Genic", "Within-LD", "Between-LD"))
trait_est.graph$Stage <- factor(trait_est.graph$Stage, levels = c("HDRW", "EYT"))
str(trait_est.graph)
head(trait_est.graph)

xbegin <- rep(c(0.62, 1, 1.62, 2, 2.62, 3, 3.62, 4), each=1000)
xend <- rep(c(0.999, 1.379, 1.999, 2.379, 2.999, 3.379, 3.999, 4.379), each=1000)
pdf("./by_chrom_sum.pdf", width=13, height=8)
graph.trait <- ggplot(trait_est.graph, aes(x=Estimates, y=Values, fill=interaction(Stage, Estimates), dodge=Stage)) +
  geom_boxplot() + geom_segment(aes(y = sum_true, yend = sum_true, x = xbegin, xend = xend), 
                                linetype="solid", size=1.2, color="#F50707") +
  geom_hline(aes(yintercept = 0), colour = "black", alpha=0.5, size=0.4, linetype="dashed") +
  scale_fill_manual(values = color_roslin)
graph.trait + optns 
dev.off()

xbegin.v <- rep(c(0.6, 1.06, 1.6, 2.06, 2.6, 3.06, 3.6, 4.06), each=1000)
xend.v <- rep(c(0.945, 1.405, 1.945, 2.405, 2.945, 3.405, 3.945, 4.405), each=1000)
pdf("./by_chrom_sum_violin.pdf", width=13, height=8)
graph.trait <- ggplot(trait_est.graph, aes(x=Estimates, y=Values, fill=interaction(Stage, Estimates), dodge=Stage)) +
  geom_violin() + geom_segment(aes(y = sum_true, yend = sum_true, x = xbegin.v, xend = xend.v), 
                                linetype="solid", size=1.2, color="#F50707") +
  geom_hline(aes(yintercept = 0), colour = "black", alpha=0.5, size=0.4, linetype="dashed") +
  scale_fill_manual(values = color_roslin)
graph.trait + optns 
dev.off()

library(coda)
HDRW.genetic.CI <- HPDinterval(as.mcmc(sum_sample.HDRW.genetic), prob=0.95)
HDRW.genic.CI <- HPDinterval(as.mcmc(sum_sample.HDRW.genic), prob=0.95)
HDRW.LD.CI <- HPDinterval(as.mcmc(sum_sample.HDRW.LD), prob=0.95)
HDRW.LD.without.CI <- HPDinterval(as.mcmc(sum_sample.HDRW.LD.without), prob=0.95)
EYT.genetic.CI <- HPDinterval(as.mcmc(sum_sample.EYT.genetic), prob=0.95)
EYT.genic.CI <- HPDinterval(as.mcmc(sum_sample.EYT.genic), prob=0.95)
EYT.LD.CI <- HPDinterval(as.mcmc(sum_sample.EYT.LD), prob=0.95)
EYT.LD.without.CI <- HPDinterval(as.mcmc(sum_sample.EYT.LD.without), prob=0.95)

(sum_true.HDRW.genetic >= HDRW.genetic.CI[1]) & (sum_true.HDRW.genetic <= HDRW.genetic.CI[2])
(sum_true.HDRW.genic >= HDRW.genic.CI[1]) & (sum_true.HDRW.genic <= HDRW.genic.CI[2])
(sum_true.HDRW.LD >= HDRW.LD.CI[1]) & (sum_true.HDRW.LD <= HDRW.LD.CI[2])
(sum_true.HDRW.LD.without >= HDRW.LD.without.CI[1]) & (sum_true.HDRW.LD.without <= HDRW.LD.without.CI[2])
(sum_true.EYT.genetic >= EYT.genetic.CI[1]) & (sum_true.EYT.genetic <= EYT.genetic.CI[2])
(sum_true.EYT.genic >= EYT.genic.CI[1]) & (sum_true.EYT.genic <= EYT.genic.CI[2])
(sum_true.EYT.LD >= EYT.LD.CI[1]) & (sum_true.EYT.LD <= EYT.LD.CI[2])
(sum_true.EYT.LD.without >= EYT.LD.without.CI[1]) & (sum_true.EYT.LD.without <= EYT.LD.without.CI[2])


#############################################
# Using chromosome values
optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=8))

color_roslin <- c("#DBAD0F", "#F2E1A5", "#C11773", "#DE99BF", "#3798D9", "#9DC4E0", "#DB760FFF", "#FAB493")

true_values <- data.frame(rep(c(HDRW.genetic_trait_true, EYT.genetic_trait_true,
                                HDRW.genic_trait_true, EYT.genic_trait_true,
                                HDRW.LD_trait_true, EYT.LD_trait_true,
                                HDRW.LD.without_trait_true, EYT.LD.without_trait_true), each=1000))
xbegin <- rep(c(0.62, 1, 1.62, 2, 2.62, 3, 3.62, 4), each=21000)
xend <- rep(c(0.999, 1.379, 1.999, 2.379, 2.999, 3.379, 3.999, 4.379), each=21000)
vector.order <- rep(rep(1:21, 8), each=1000)
values.temp <- cbind(true_values, xbegin, xend, vector.order)
colnames(values.temp) <- c("true.value", "xbegin", "xend", "vector.order")
values.order <- values.temp[order(values.temp$vector.order, values.temp$xbegin),]
head(values.order)

all.chr <- as.data.frame(t(cbind(sample.HDRW.genetic, sample.HDRW.genic,
                                 sample.HDRW.LD, sample.HDRW.LD.without,
                                 sample.EYT.genetic, sample.EYT.genic, 
                                 sample.EYT.LD, sample.EYT.LD.without)))
by.chr <- as.data.frame(unlist(all.chr))
colnames(by.chr) <- c("Values")
by.chr$Chr <- factor(rep(1:21, each=8000))
by.chr$Stage <- rep(rep(c("HDRW", "EYT"), each=4000), 21)
by.chr$Estimates <- rep(rep(c("Genetic", "Genic", "Within-LD", "Between-LD"), each=1000), 42)
by.chr$Estimates <- factor(by.chr$Estimates, levels = c("Genetic", "Genic", "Within-LD", "Between-LD"))
by.chr$Stage <- factor(by.chr$Stage, levels = c("HDRW", "EYT"))
str(by.chr)
head(by.chr)

pdf("./graph_byChr.pdf", width=28, height=30)
graph.trait <- ggplot(by.chr, aes(x=Estimates, y=Values, fill=interaction(Stage, Estimates), dodge=Stage)) +
  geom_boxplot() +
  geom_segment(aes(y = values.order$true.value, yend = values.order$true.value, 
                   x = values.order$xbegin, xend = values.order$xend), 
                   linetype="solid", size=1.2, color="#F50707") +
  facet_wrap( ~ Chr, ncol = 4) +
  geom_hline(aes(yintercept = 0), colour = "black", alpha=0.5, size=0.4, linetype="dashed") +
  scale_fill_manual(values = color_roslin)
graph.trait + optns 
dev.off()

true_values <- data.frame(rep(c(HDRW.genetic_trait_true, EYT.genetic_trait_true,
                                HDRW.genic_trait_true, EYT.genic_trait_true,
                                HDRW.LD_trait_true, EYT.LD_trait_true,
                                HDRW.LD.without_trait_true, EYT.LD.without_trait_true), each=1000))
xbegin.v <- rep(c(0.6, 1.06, 1.6, 2.06, 2.6, 3.06, 3.6, 4.06), each=21000)
xend.v <- rep(c(0.945, 1.405, 1.945, 2.405, 2.945, 3.405, 3.945, 4.405), each=21000)
vector.order <- rep(rep(1:21, 8), each=1000)
values.temp <- cbind(true_values, xbegin.v, xend.v, vector.order)
colnames(values.temp) <- c("true.value", "xbegin.v", "xend.v", "vector.order")
values.order <- values.temp[order(values.temp$vector.order, values.temp$xbegin.v),]
head(values.order)
pdf("./graph_byChr_violin.pdf", width=28, height=30)
graph.trait <- ggplot(by.chr, aes(x=Estimates, y=Values, fill=interaction(Stage, Estimates), dodge=Stage)) +
  geom_violin() +
  geom_segment(aes(y = values.order$true.value, yend = values.order$true.value, 
                   x = values.order$xbegin.v, xend = values.order$xend.v), 
               linetype="solid", size=1.2, color="#F50707") +
  facet_wrap( ~ Chr, ncol = 4) +
  geom_hline(aes(yintercept = 0), colour = "black", alpha=0.5, size=0.4, linetype="dashed") +
  scale_fill_manual(values = color_roslin)
graph.trait + optns 
dev.off()

library(coda)
HDRW.genetic.CI <- apply(sample.HDRW.genetic, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
HDRW.genic.CI <- apply(sample.HDRW.genic, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
HDRW.LD.CI <- apply(sample.HDRW.LD, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
HDRW.LD.without.CI <- apply(sample.HDRW.LD.without, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
EYT.genetic.CI <- apply(sample.EYT.genetic, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
EYT.genic.CI <- apply(sample.EYT.genic, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
EYT.LD.CI <- apply(sample.EYT.LD, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))
EYT.LD.without.CI <- apply(sample.EYT.LD.without, 1, function(x) HPDinterval(as.mcmc(x), prob=0.95))

(HDRW.genetic_trait_true >= HDRW.genetic.CI[1,]) & (HDRW.genetic_trait_true <= HDRW.genetic.CI[2,])
(HDRW.genic_trait_true >= HDRW.genic.CI[1,]) & (HDRW.genic_trait_true <= HDRW.genic.CI[2,])
(HDRW.LD_trait_true >= HDRW.LD.CI[1,]) & (HDRW.LD_trait_true <= HDRW.LD.CI[2,])
(HDRW.LD.without_trait_true >= HDRW.LD.without.CI[1,]) & (HDRW.LD.without_trait_true <= HDRW.LD.without.CI[2,])
(EYT.genetic_trait_true >= EYT.genetic.CI[1,]) & (EYT.genetic_trait_true <= EYT.genetic.CI[2,])
(EYT.genic_trait_true >= EYT.genic.CI[1,]) & (EYT.genic_trait_true <= EYT.genic.CI[2,])
(EYT.LD_trait_true >= EYT.LD.CI[1,]) & (EYT.LD_trait_true <= EYT.LD.CI[2,])
(EYT.LD.without_trait_true >= EYT.LD.without.CI[1,]) & (EYT.LD.without_trait_true <= EYT.LD.without.CI[2,])


# "#A9A9A9" -> silver

##############################################################################################################
##############################################################################################################

#############################################
#### Graph - Genome ####
library(ggplot2)
rm(list=ls())

load("byChrom_genome_QTLs_EYT.RData")
EYT.genetic_genome_QTLs <- unlist(chr.genetic)
EYT.genic_genome_QTLs <- unlist(chr.genic)
EYT.LD_genome_QTLs <- unlist(chr.LD)
EYT.LD.without_genome_QTLs <- unlist(chr.LD.without)

load("byChrom_genome_QTLs_HDRW.RData")
HDRW.genetic_genome_QTLs <- unlist(chr.genetic)
HDRW.genic_genome_QTLs <- unlist(chr.genic)
HDRW.LD_genome_QTLs <- unlist(chr.LD)
HDRW.LD.without_genome_QTLs <- unlist(chr.LD.without)

load("byChrom_genome_SNPs_EYT.RData")
EYT.genetic_genome_SPNs <- unlist(chr.genetic)
EYT.genic_genome_SNPs <- unlist(chr.genic)
EYT.LD_genome_SNPs <- unlist(chr.LD)
EYT.LD.without_genome_SNPs <- unlist(chr.LD.without)

load("byChrom_genome_SNPs_HDRW.RData")
HDRW.genetic_genome_SNPs <- unlist(chr.genetic)
HDRW.genic_genome_SNPs <- unlist(chr.genic)
HDRW.LD_genome_SNPs <- unlist(chr.LD)
HDRW.LD.without_genome_SNPs <- unlist(chr.LD.without)




