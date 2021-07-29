
rm(list = ls())
library(AlphaSimR)             # Version 1.0.1

QTLs = c(10, 100, 1000)
nSNPs = 500
scenGxY = c("noGxY", "GxY")

for(r in 1:10){
  for(intGxY in scenGxY){
    for(nQTLs in QTLs){
      source("BURNIN.R") 
      source("BASE.R")
    }
  }
}





