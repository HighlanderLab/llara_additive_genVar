
library(AlphaSimR)
library(tidyverse)
library(dplyr)

#### Full Bayes ####
for(simulation in 1:10){
  REP = paste0("_Sim", simulation)
  source("Estimates_full.R") 
}

#### Empirical Bayes ####
for(simulation in 1:10){
  REP = paste0("_Sim", simulation)
  source("Estimates_empirical.R") 
}


