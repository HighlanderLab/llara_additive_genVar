
#######################################################################
#            Creating Graphs Empirical x Full Bayesian                #
#                 Simulated wheat breeding program                    #
#                   Leticia Lara - February, 2020                     #
#######################################################################

library(ggplot2)
library(coda)
library(openxlsx)
library(tidyr)
library(ggridges)
library(GGally)
library(ggpubr)
library(AlphaSimR)           # Version 0.11.0


rm(list=ls())
####################
### Sumarizing results - Table
data_emp.graph1 <- NULL
data_full.graph1 <- NULL
data_emp.graph2 <- list(year16 = list(), year17 = list(), year18 = list(),
                        year19 = list(), year20 = list(), year21 = list()) 
data_full.graph2 <- list(year16 = list(), year17 = list(), year18 = list(),
                         year19 = list(), year20 = list(), year21 = list()) 
i=1

# Empirical Bayes
for(year in 16:21){
  load(paste0("./empiricalBayes/1000samples.empiricalBayes_year",year,".RData"))

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
  
  table_results <- data.frame(Stage = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"),
                              True.genetic.var = c(stage.year.Parents$var, stage.year.F1$var, stage.year.HDRW$var,
                                                   stage.year.PYT$var, stage.year.AYT$var, EYT1.true.Genetic),
                              Estimated.genetic.var = c(mean(Parents.Genetic.var), mean(F1.Genetic.var),
                                                        mean(HDRW.Genetic.var), mean(PYT.Genetic.var),
                                                        mean(AYT.Genetic.var), mean(EYT1.Genetic.var)),
                              Est.genetic.sd.var = c(sd(Parents.Genetic.var), sd(F1.Genetic.var),
                                                     sd(HDRW.Genetic.var), sd(PYT.Genetic.var),
                                                     sd(AYT.Genetic.var), sd(EYT1.Genetic.var)),
                              Lower.genetic = Genetic.CI[,1],
                              Upper.genetic = Genetic.CI[,2],
                              True.genic.var = c(stage.year.Parents$genicG, stage.year.F1$genicG, stage.year.HDRW$genicG,
                                                 stage.year.PYT$genicG, stage.year.AYT$genicG, EYT1.true.Genic),
                              Estimated.genic.var = c(mean(Parents.Genic.var), mean(F1.Genic.var),
                                                      mean(HDRW.Genic.var), mean(PYT.Genic.var),
                                                      mean(AYT.Genic.var), mean(EYT1.Genic.var)),
                              Est.genic.sd.var = c(sd(Parents.Genic.var), sd(F1.Genic.var),
                                                   sd(HDRW.Genic.var), sd(PYT.Genic.var),
                                                   sd(AYT.Genic.var), sd(EYT1.Genic.var)),
                              Lower.genic = Genic.CI[,1],
                              Upper.genic = Genic.CI[,2])
  data_emp.graph1 <- rbind(data_emp.graph1, table_results)
  data_emp.graph2[[i]] <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                                     PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                                     Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                                   PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                                     Stage = c(rep("Parents",1000), rep("F1",1000), rep("HDRW",1000),
                                               rep("PYT",1000), rep("AYT",1000), rep("EYT",1000)))
  data_emp.graph2[[i]]$Stage <- factor(data_emp.graph2[[i]]$Stage, levels = c("EYT", "AYT", "PYT", "HDRW", "F1", "Parents"))
  i=i+1
}

# Full Bayes
for(year in 16:21){
  load(paste0("./fullBayes/1000samples.fullBayes_year",year,".RData"))

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
  
  table_results <- data.frame(Stage = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"),
                              True.genetic.var = c(stage.year.Parents$var, stage.year.F1$var, stage.year.HDRW$var,
                                                   stage.year.PYT$var, stage.year.AYT$var, EYT1.true.Genetic),
                              Estimated.genetic.var = c(mean(Parents.Genetic.var), mean(F1.Genetic.var),
                                                        mean(HDRW.Genetic.var), mean(PYT.Genetic.var),
                                                        mean(AYT.Genetic.var), mean(EYT1.Genetic.var)),
                              Est.genetic.sd.var = c(sd(Parents.Genetic.var), sd(F1.Genetic.var),
                                                     sd(HDRW.Genetic.var), sd(PYT.Genetic.var),
                                                     sd(AYT.Genetic.var), sd(EYT1.Genetic.var)),
                              Lower.genetic = Genetic.CI[,1],
                              Upper.genetic = Genetic.CI[,2],
                              True.genic.var = c(stage.year.Parents$genicG, stage.year.F1$genicG, stage.year.HDRW$genicG,
                                                 stage.year.PYT$genicG, stage.year.AYT$genicG, EYT1.true.Genic),
                              Estimated.genic.var = c(mean(Parents.Genic.var), mean(F1.Genic.var),
                                                      mean(HDRW.Genic.var), mean(PYT.Genic.var),
                                                      mean(AYT.Genic.var), mean(EYT1.Genic.var)),
                              Est.genic.sd.var = c(sd(Parents.Genic.var), sd(F1.Genic.var),
                                                   sd(HDRW.Genic.var), sd(PYT.Genic.var),
                                                   sd(AYT.Genic.var), sd(EYT1.Genic.var)),
                              Lower.genic = Genic.CI[,1],
                              Upper.genic = Genic.CI[,2])
  data_full.graph1 <- rbind(data_full.graph1, table_results)
  data_full.graph2[[i]] <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                                      PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                                     Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                                   PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                                     Stage = c(rep("Parents",1000), rep("F1",1000), rep("HDRW",1000),
                                               rep("PYT",1000), rep("AYT",1000), rep("EYT",1000)))
  data_full.graph2[[i]]$Stage <- factor(data_full.graph2[[i]]$Stage, levels = c("EYT", "AYT", "PYT", "HDRW", "F1", "Parents"))
  i=i+1
}


####################
# Graph: Empirical x Full x True
True.value <- c(data_emp.graph1[,2], data_emp.graph1[,7])
Empirical <- c(data_emp.graph1[,3], data_emp.graph1[,8])
Full <- c(data_full.graph1[,3], data_full.graph1[,8])
Variance <- c(rep("Genetic", 36), rep("Genic", 36))
Stage <- rep(data_emp.graph1[,1], 2)
Year <- rep(c(rep(16:21,each=6)),2)

temp.data1 <- data.frame(Stage, Year, Variance)
temp.data1$Stage <- factor(temp.data1$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
data.graph1 <- data.frame(True.value, Empirical, Full)
head(data.graph1)

optns <- theme_bw(base_size = 18, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")

data_genetic.graph1 <- data.graph1[1:36,]
data_genic.graph1 <- data.graph1[37:72,]
stage.data1 <- temp.data1[1:36,1:2]

#Genetic
pdf("./graph_Genetic.pdf", width=11, height=10)
graph1.genetic <- ggpairs(data_genetic.graph1, aes(colour = stage.data1$Stage))

for(i in 1:graph1.genetic$nrow) {
  for(j in 1:graph1.genetic$ncol){
    graph1.genetic[i,j] <- graph1.genetic[i,j] + 
      scale_fill_manual(values = color_ci) +
      scale_color_manual(values = color_rgb) 
  }
}
graph1.genetic + optns + ggtitle("Genetic variance")
dev.off()

#Genic
pdf("./graph_Genic.pdf", width=11, height=10)
graph1.genic <- ggpairs(data_genic.graph1, aes(colour = stage.data1$Stage))

for(i in 1:graph1.genic$nrow) {
  for(j in 1:graph1.genic$ncol){
    graph1.genic[i,j] <- graph1.genic[i,j] + 
      scale_fill_manual(values = color_ci) +
      scale_color_manual(values = color_rgb)  
  }
}
graph1.genic + optns + ggtitle("Genic variance") 
dev.off()



##############################################################################################################
##############################################################################################################

# Graph: Empirical x Full
Empirical <- c(data_emp.graph1[,3], data_emp.graph1[,8])
Full <- c(data_full.graph1[,3], data_full.graph1[,8])
Variance <- c(rep("Genetic", 36), rep("Genic", 36))
Stage <- rep(data_emp.graph1[,1], 2)
Year <- rep(c(rep(16:21,each=6)),2)

temp.data1 <- data.frame(Stage, Year, Variance)
temp.data1$Stage <- factor(temp.data1$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
data.graph1 <- data.frame(Empirical, Full)
head(data.graph1)

optns <- theme_bw(base_size = 18, base_family = "sans") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")

data_genetic.graph1 <- data.graph1[1:36,]
data_genic.graph1 <- data.graph1[37:72,]
stage.data1 <- temp.data1[1:36,1:2]

#Genetic
pdf("./graph_Genetic.pdf", width=11, height=10)
graph1.genetic <- ggpairs(data_genetic.graph1, aes(colour = stage.data1$Stage))

for(i in 1:graph1.genetic$nrow) {
  for(j in 1:(graph1.genetic$ncol)){
    graph1.genetic[i,j] <- graph1.genetic[i,j] + 
      scale_fill_manual(values = color_ci) +
      scale_color_manual(values = color_rgb) 
  }
}
graph1.genetic + optns + ggtitle("Genetic variance") 
dev.off()

#Genic
pdf("./graph_Genic.pdf", width=11, height=10)
graph1.genic <- ggpairs(data_genic.graph1, aes(colour = stage.data1$Stage))

for(i in 1:graph1.genic$nrow) {
  for(j in 1:graph1.genic$ncol){
    graph1.genic[i,j] <- graph1.genic[i,j] + 
      scale_fill_manual(values = color_ci) +
      scale_color_manual(values = color_rgb)  
  }
}
graph1.genic + optns + ggtitle("Genic variance") 
dev.off()


##############################################################################################################
##############################################################################################################

####################
# Graph: Correlations Empirical x Full for Genetic, Genic, Mean, and SD

optns <- theme_bw(base_size = 18, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")
shape.opt <- c("15", "16", "17", "23", "8", "11")

# Genetic - Mean
head(data_emp.graph1)
Year <- rep(16:21,each=6)
data.graph2.1 <- data.frame(data_emp.graph1$Stage, Year, data_emp.graph1$Estimated.genetic.var, 
                            data_full.graph1$Estimated.genetic.var)
colnames(data.graph2.1) <- c("Stage", "Year", "Value.emp", "Value.full")
data.graph2.1$Year <- as.factor(data.graph2.1$Year)
data.graph2.1$Stage <- factor(data.graph2.1$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
head(data.graph2.1)
graph2.1 <- ggplot(data.graph2.1, aes(Value.full, Value.emp)) + 
  geom_point(size=2) + coord_cartesian(ylim = c(0.05, 0.10), xlim = c(0.05, 0.10)) + 
  geom_abline(aes(intercept = 0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") +
  labs(x = "Full Bayes", y = "Empirical Bayes") + ggtitle("Posterior Mean - Genetic variance") + optns
  
# Genetic - SD
data.graph2.2 <- data.frame(data_emp.graph1$Stage, Year, data_emp.graph1$Est.genetic.sd.var, 
                            data_full.graph1$Est.genetic.sd.var)
colnames(data.graph2.2) <- c("Stage", "Year", "Value.emp", "Value.full")
data.graph2.2$Year <- as.factor(data.graph2.2$Year)
data.graph2.2$Stage <- factor(data.graph2.2$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
head(data.graph2.2)
graph2.2 <- ggplot(data.graph2.2, aes(Value.full, Value.emp)) + 
  geom_point(size=2) + coord_cartesian(ylim = c(0.003, 0.03), xlim = c(0.003, 0.03)) + 
  geom_abline(aes(intercept = 0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") + 
  labs(x = "Full Bayes", y = "Empirical Bayes") + ggtitle("Posterior Standard Deviation - Genetic variance") + optns

# Genic - Mean
data.graph2.3 <- data.frame(data_emp.graph1$Stage, Year, data_emp.graph1$Estimated.genic.var, 
                            data_full.graph1$Estimated.genic.var)
colnames(data.graph2.3) <- c("Stage", "Year", "Value.emp", "Value.full")
data.graph2.3$Year <- as.factor(data.graph2.3$Year)
data.graph2.3$Stage <- factor(data.graph2.3$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
head(data.graph2.3)
graph2.3 <- ggplot(data.graph2.3, aes(Value.full, Value.emp)) + 
  geom_point(size=2) + coord_cartesian(ylim = c(0.04, 0.10), xlim = c(0.04, 0.10)) + 
  geom_abline(aes(intercept = 0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") + 
  labs(x = "Full Bayes", y = "Empirical Bayes") + ggtitle("Posterior Mean - Genic variance") + optns

# Genic - SD
data.graph2.4 <- data.frame(data_emp.graph1$Stage, Year, data_emp.graph1$Est.genic.sd.var, 
                            data_full.graph1$Est.genic.sd.var)
colnames(data.graph2.4) <- c("Stage", "Year", "Value.emp", "Value.full")
data.graph2.4$Year <- as.factor(data.graph2.4$Year)
data.graph2.4$Stage <- factor(data.graph2.4$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
head(data.graph2.4)
graph2.4 <- ggplot(data.graph2.4, aes(Value.full, Value.emp)) + 
  geom_point(size=2) + coord_cartesian(ylim = c(0, 0.02), xlim = c(0, 0.02)) + 
  geom_abline(aes(intercept = 0, slope=1), colour = "red", alpha=0.5, size=0.4, linetype="dashed") + 
  labs(x = "Full Bayes", y = "Empirical Bayes") + ggtitle("Posterior Standard Deviation - Genic variance") + optns


# All
pdf("./empirical_X_full.pdf", width=16, height=10)
graph2 <- ggarrange(graph2.1, graph2.2, graph2.3, graph2.4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
graph2
dev.off()


##############################################################################################################
##############################################################################################################

