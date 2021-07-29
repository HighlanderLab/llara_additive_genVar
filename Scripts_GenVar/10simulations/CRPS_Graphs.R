
#######################################################################
#          CRPS analysis and results from 10 simulations              #
#                 Simulated wheat breeding program                    #
#                    Leticia Lara - June, 2020                        #
#######################################################################

library(verification)
library(scoringRules)
library(AlphaSimR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)

rm(list=ls())


##########################################
##### Sumarizing results - FullBayes #####
StageName <- c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT")
results.summary_Full <- results.sample_Full <- NULL
for(sim in 1:10){
  load(paste0("./AlphaSimR/Results/Burnin_Sim", sim, ".RData"))
  load(paste0("./AlphaSimR/Results/all_stages_and_years_Sim", sim, ".RData"))
  
  for(year in 16:21){
    load(paste0("./AlphaBayes-Estimation/Simulation", sim, "/fullBayes/Estimates.fullBayes_year", year, "_Sim", sim,".RData"))
    
    Genetic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(F1.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(HDRW.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(PYT.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(AYT.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(EYT1.Genetic.var), prob=0.95))
    Genic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(F1.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(HDRW.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(PYT.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(AYT.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(EYT1.Genic.var), prob=0.95))
    
    table_results <- data.frame(Stage = StageName,
                                Year = rep(year, length(StageName)),
                                Sim = rep(sim, length(StageName)),
                                True.genetic.var = c(stage.year.Parents$var, stage.year.F1$var, stage.year.HDRW$var,
                                                     stage.year.PYT$var, stage.year.AYT$var, EYT1.true.Genetic),
                                Estimated.genetic.var = c(mean(Parents.Genetic.var), mean(F1.Genetic.var),
                                                          mean(HDRW.Genetic.var), mean(PYT.Genetic.var),
                                                          mean(AYT.Genetic.var), mean(EYT1.Genetic.var)),
                                SD.genetic.var = c(sd(Parents.Genetic.var), sd(F1.Genetic.var), sd(HDRW.Genetic.var), 
                                                          sd(PYT.Genetic.var), sd(AYT.Genetic.var), sd(EYT1.Genetic.var)),
                                Lower.genetic = Genetic.CI[,1],
                                Upper.genetic = Genetic.CI[,2],
                                True.genic.var = c(stage.year.Parents$genicG, stage.year.F1$genicG, stage.year.HDRW$genicG,
                                                   stage.year.PYT$genicG, stage.year.AYT$genicG, EYT1.true.Genic),
                                Estimated.genic.var = c(mean(Parents.Genic.var), mean(F1.Genic.var), mean(HDRW.Genic.var), 
                                                        mean(PYT.Genic.var), mean(AYT.Genic.var), mean(EYT1.Genic.var)),
                                SD.genic.var = c(sd(Parents.Genic.var), sd(F1.Genic.var), sd(HDRW.Genic.var), 
                                                 sd(PYT.Genic.var), sd(AYT.Genic.var), sd(EYT1.Genic.var)),
                                Lower.genic = Genic.CI[,1],
                                Upper.genic = Genic.CI[,2])
    results.summary_Full <- rbind(results.summary_Full, table_results)
    
    sample <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                         PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                         Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                       PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                         Stage = rep(StageName, each = 900),
                         Year = year,
                         Sim =  sim)
    results.sample_Full <- rbind(results.sample_Full, sample)
    
    Parents = parents.c[rep(16:21, each=70) == year]
    F1 = F1c[rep(16:21, each=100) == year]
    HDRW = HDRWc[c(rep(16:20, each=10000), rep(21, 500)) == year]
    PYT = PYTc[rep(16:21, each=500) == year]
    AYT = AYTc[rep(16:21, each=50) == year]
    EYT.year <- EYTc[rep(16:21, each=20) == year]
    EYT = EYT.year[1:10]
    pheno.Stage <- list(Parents = data.frame(Parents@pheno),
                        F1 = data.frame(F1@pheno),
                        HDRW = data.frame(HDRW@pheno),
                        PYT = data.frame(PYT@pheno),
                        AYT = data.frame(AYT@pheno),
                        EYT = data.frame(EYT@pheno))
    if((sim==1)&(year==16)){
      pheno_Full = pheno.Stage
    } else{
      pheno_Full <- mapply(cbind, pheno_Full, pheno.Stage, SIMPLIFY = FALSE) 
    }
  }
}
col_names <- c(paste0("Sim", rep(1:10, each=6), ".year", c(16:21)))
pheno_Full <- lapply(pheno_Full, setNames, nm = col_names)


###############################################
##### Sumarizing results - EmpiricalBayes #####
results.summary_Empirical <- results.sample_Empirical <- NULL
for(sim in 1:10){
  load(paste0("./AlphaSimR/Results/Burnin_Sim", sim, ".RData"))
  load(paste0("./AlphaSimR/Results/all_stages_and_years_Sim", sim, ".RData"))
  
  for(year in 16:21){
    load(paste0("./AlphaBayes-Estimation/Simulation", sim, "/empiricalBayes/Estimates.empiricalBayes_year", year, "_Sim", sim, ".RData"))
    
    Genetic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(F1.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(HDRW.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(PYT.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(AYT.Genetic.var), prob=0.95), 
                        HPDinterval(as.mcmc(EYT1.Genetic.var), prob=0.95))
    Genic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(F1.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(HDRW.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(PYT.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(AYT.Genic.var), prob=0.95), 
                      HPDinterval(as.mcmc(EYT1.Genic.var), prob=0.95))
    
    table_results <- data.frame(Stage = StageName,
                                Year = rep(year, length(StageName)),
                                Sim = rep(sim, length(StageName)),
                                True.genetic.var = c(stage.year.Parents$var, stage.year.F1$var, stage.year.HDRW$var,
                                                     stage.year.PYT$var, stage.year.AYT$var, EYT1.true.Genetic),
                                Estimated.genetic.var = c(mean(Parents.Genetic.var), mean(F1.Genetic.var),
                                                          mean(HDRW.Genetic.var), mean(PYT.Genetic.var),
                                                          mean(AYT.Genetic.var), mean(EYT1.Genetic.var)),
                                SD.genetic.var = c(sd(Parents.Genetic.var), sd(F1.Genetic.var), sd(HDRW.Genetic.var), 
                                                   sd(PYT.Genetic.var), sd(AYT.Genetic.var), sd(EYT1.Genetic.var)),
                                Lower.genetic = Genetic.CI[,1],
                                Upper.genetic = Genetic.CI[,2],
                                True.genic.var = c(stage.year.Parents$genicG, stage.year.F1$genicG, stage.year.HDRW$genicG,
                                                   stage.year.PYT$genicG, stage.year.AYT$genicG, EYT1.true.Genic),
                                Estimated.genic.var = c(mean(Parents.Genic.var), mean(F1.Genic.var),
                                                        mean(HDRW.Genic.var), mean(PYT.Genic.var),
                                                        mean(AYT.Genic.var), mean(EYT1.Genic.var)),
                                SD.genic.var = c(sd(Parents.Genic.var), sd(F1.Genic.var), sd(HDRW.Genic.var), 
                                                 sd(PYT.Genic.var), sd(AYT.Genic.var), sd(EYT1.Genic.var)),
                                Lower.genic = Genic.CI[,1],
                                Upper.genic = Genic.CI[,2])
    results.summary_Empirical <- rbind(results.summary_Empirical, table_results)

    sample <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                         PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                         Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                       PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                         Stage = rep(StageName, each = 900),
                         Year = year,
                         Sim =  sim)
    results.sample_Empirical <- rbind(results.sample_Empirical, sample)
    
    Parents = parents.c[rep(16:21, each=70) == year]
    F1 = F1c[rep(16:21, each=100) == year]
    HDRW = HDRWc[c(rep(16:20, each=10000), rep(21, 500)) == year]
    PYT = PYTc[rep(16:21, each=500) == year]
    AYT = AYTc[rep(16:21, each=50) == year]
    EYT.year <- EYTc[rep(16:21, each=20) == year]
    EYT = EYT.year[1:10]
    pheno.Stage <- list(Parents = data.frame(Parents@pheno),
                        F1 = data.frame(F1@pheno),
                        HDRW = data.frame(HDRW@pheno),
                        PYT = data.frame(PYT@pheno),
                        AYT = data.frame(AYT@pheno),
                        EYT = data.frame(EYT@pheno))
    if((sim==1)&(year==16)){
      pheno_Empirical = pheno.Stage
    } else{
      pheno_Empirical <- mapply(cbind, pheno_Empirical, pheno.Stage, SIMPLIFY = FALSE) 
    }
  }
}
pheno_Empirical <- lapply(pheno_Empirical, setNames, nm = col_names)
#save(results.summary_Full, results.sample_Full, results.summary_Empirical, 
#     results.sample_Empirical, pheno_Full, pheno_Empirical, file = "results_Simulations.RData")


########################################
##### CRPS analysis - scoringRules #####
rm(list=ls())
load("results_Simulations.RData")

StageName <- c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT")
crps.genic_Full <- crps.genic_Empirical <- NULL
crps.genetic_Full <- crps.genetic_Empirical <- NULL
for(sim in 1:10){
  for(year in 16:21){
    for(stage in 1:length(StageName)){
      TrueGenic.full <- results.summary_Full %>% 
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$True.genic.var
      TrueGenetic.full <- results.summary_Full %>% 
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$True.genetic.var

      sampleGenic.full <- results.sample_Full %>%
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$Genic.var
      sampleGenetic.full <- results.sample_Full %>%
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$Genetic.var

      crpsGenic.full <- crps_sample(TrueGenic.full, dat = sampleGenic.full)
      crpsGenetic.full <- crps_sample(TrueGenetic.full, dat = sampleGenetic.full)
      
      TrueGenic.emp <- results.summary_Empirical %>% 
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$True.genic.var
      TrueGenetic.emp <- results.summary_Empirical %>% 
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$True.genetic.var
      
      sampleGenic.emp <- results.sample_Empirical %>%
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$Genic.var
      sampleGenetic.emp <- results.sample_Empirical %>%
        filter((Stage == StageName[stage])&(Year == year)&(Sim == sim)) %>% .$Genetic.var
      
      crpsGenic.emp <- crps_sample(TrueGenic.emp, dat = sampleGenic.emp)
      crpsGenetic.emp <- crps_sample(TrueGenetic.emp, dat = sampleGenetic.emp)
      
      crps.genic_Full <- c(crps.genic_Full, crpsGenic.full)
      crps.genic_Empirical <- c(crps.genic_Empirical, crpsGenic.emp)
      crps.genetic_Full <- c(crps.genetic_Full, crpsGenetic.full)
      crps.genetic_Empirical <- c(crps.genetic_Empirical, crpsGenetic.emp)
    }
  }
}
crps.scoringRules <- cbind(results.summary_Full[,c(1:3)], crps.genic_Full, crps.genic_Empirical, 
                           crps.genetic_Full, crps.genetic_Empirical)
head(crps.scoringRules)

tapply(crps.scoringRules$crps.genic_Full, crps.scoringRules$Stage, mean)
tapply(crps.scoringRules$crps.genic_Empirical, crps.scoringRules$Stage, mean)
tapply(crps.scoringRules$crps.genetic_Full, crps.scoringRules$Stage, mean)
tapply(crps.scoringRules$crps.genetic_Empirical, crps.scoringRules$Stage, mean)

tapply(crps.scoringRules$crps.genic_Full, crps.scoringRules$Stage, sd)
tapply(crps.scoringRules$crps.genic_Empirical, crps.scoringRules$Stage, sd)
tapply(crps.scoringRules$crps.genetic_Full, crps.scoringRules$Stage, sd)
tapply(crps.scoringRules$crps.genetic_Empirical, crps.scoringRules$Stage, sd)

# Graph
for(sim in 1:10){
  crps.sim <- filter(crps.scoringRules, Sim == sim)
  crps.data <- gather(crps.sim, Estimates, Value, crps.genic_Full:crps.genetic_Empirical, factor_key = TRUE)
  crps.data$Analysis <- as.factor(unlist(lapply(strsplit(as.character(crps.data$Estimates), "_"), "[", 2)))
  crps.data$Estimates <- as.factor(unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(crps.data$Estimates), 
                                                                                 "_"), "[", 1)), "[.]"), "[", 2)))
  crps.data$Stage <- factor(crps.data$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
  
  graph.sim <- ggplot(crps.data, aes(y = Value, x = Year, group = interaction(Estimates, Analysis))) +
    geom_line(aes(color = Analysis, linetype = Estimates)) + xlab("Stage") + ylab("CRPS") + 
    facet_wrap(~ Stage) + coord_cartesian(ylim = c(0, 0.058)) +
    scale_color_manual(values = c("#9d63eb", "#78cc8a")) + ggtitle(paste0("Simulation ", sim)) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  update_geom_defaults("point", list(size = 1))
  theme_set(theme_grey(base_size = 6))
  ggsave(filename = paste0("myGraph", sim, ".png"), plot=graph.sim)
}
# convert -delay 50 myGraph*.png CRPS_simulations.gif


########################################
##### CRPS analysis - verification #####
rm(list=ls())
load("results_Simulations.RData")

StageName <- c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT")
crps.genic_Full <- crps.genic_Empirical <- NULL
crps.genetic_Full <- crps.genetic_Empirical <- NULL
for(sim in 1:10){
  for(year in 16:21){
    for(stage in 1:length(StageName)){
      Est.full <- results.summary_Full %>% filter((Stage == StageName[stage])&(Year == year)&(Sim == sim))
      crpsGenic.full <- verification::crps(Est.full$True.genic.var, c(Est.full$Estimated.genic.var, Est.full$SD.genic.var))
      crpsGenetic.full <- verification::crps(Est.full$True.genetic.var, c(Est.full$Estimated.genetic.var, Est.full$SD.genetic.var))
      
      Est.emp <- results.summary_Empirical %>% filter((Stage == StageName[stage])&(Year == year)&(Sim == sim))
      Obs.emp <- as.matrix(pheno_Empirical[[stage]][,paste0("Sim", sim, ".year", year)])
      crpsGenic.emp <- verification::crps(Est.emp$True.genic.var, c(Est.emp$Estimated.genic.var, Est.emp$SD.genic.var))
      crpsGenetic.emp <- verification::crps(Est.emp$True.genetic.var, c(Est.emp$Estimated.genetic.var, Est.emp$SD.genetic.var))
      
      crps.genic_Full <- c(crps.genic_Full, crpsGenic.full$CRPS)
      crps.genic_Empirical <- c(crps.genic_Empirical, crpsGenic.emp$CRPS)
      crps.genetic_Full <- c(crps.genetic_Full, crpsGenetic.full$CRPS)
      crps.genetic_Empirical <- c(crps.genetic_Empirical, crpsGenetic.emp$CRPS)
    }
  }
}
crps.verification <- cbind(results.summary_Full[,c(1:3)], crps.genic_Full, crps.genic_Empirical, 
                           crps.genetic_Full, crps.genetic_Empirical)
head(crps.verification)

# Graph
optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size = 1))
theme_set(theme_grey(base_size = 6))

# Correlation graph
graph.genic <- ggplot(crps.verification, aes(crps.genic_Full, crps.genic_Empirical)) + 
  geom_point(size = 2) + labs(x = "Full Bayesian", y = "Empirical Bayesian") + 
  geom_abline(aes(intercept=0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") +
  ggtitle("Genic variance") + optns
graph.genetic <- ggplot(crps.verification, aes(crps.genetic_Full, crps.genetic_Empirical)) + 
  geom_point(size = 2) + labs(x = "Full Bayesian", y = "Empirical Bayesian") + 
  geom_abline(aes(intercept=0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") +
  ggtitle("Genetic variance") + optns

cor(crps.verification$crps.genic_Full, crps.genic_Empirical)
cor(crps.verification$crps.genetic_Full, crps$crps.genetic_Empirical)

for(sim in 1:10){
  crps.sim <- filter(crps.verification, Sim == sim)
  crps.data <- gather(crps.sim, Estimates, Value, crps.genic_Full:crps.genetic_Empirical, factor_key = TRUE)
  crps.data$Analysis <- as.factor(unlist(lapply(strsplit(as.character(crps.data$Estimates), "_"), "[", 2)))
  crps.data$Estimates <- as.factor(unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(crps.data$Estimates), 
                                                                                 "_"), "[", 1)), "[.]"), "[", 2)))
  crps.data$Stage <- factor(crps.data$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
  
  graph.sim <- ggplot(crps.data, aes(y = Value, x = Year, group = interaction(Estimates, Analysis))) +
    geom_line(aes(color = Analysis, linetype = Estimates)) + xlab("Stage") + ylab("CRPS") + 
    facet_wrap(~ Stage) + 
    scale_color_manual(values = c("#9d63eb", "#78cc8a")) + ggtitle(paste0("Simulation ", sim)) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  update_geom_defaults("point", list(size = 1))
  theme_set(theme_grey(base_size = 6))
  ggsave(filename = paste0("myGraph_verification", sim, ".png"), plot = graph.sim)
}
# convert -delay 50 myGraph_verification*.png CRPS_simulationsVerification.gif

cor(crps.verification$crps.genetic_Full, crps.scoringRules$crps.genetic_Full)
cor(crps.verification$crps.genetic_Empirical, crps.scoringRules$crps.genetic_Empirical)
cor(crps.verification$crps.genic_Full, crps.scoringRules$crps.genic_Full)
cor(crps.verification$crps.genic_Empirical, crps.scoringRules$crps.genic_Empirical)



##############################################################################################################
##############################################################################################################

