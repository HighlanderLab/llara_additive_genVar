
#######################################################################
#               Creating Graphs Principal Components                  #
#                 Simulated wheat breeding program                    #
#                    Leticia Lara - March, 2020                       #
#######################################################################

library(ggplot2)
library(coda)
library(openxlsx)
library(tidyr)
library(ggridges)
library(GGally)


rm(list=ls())
####################
### Sumarizing results - Table
pc <- c(10, 50, 100, 500, 1000, 2000, 3410)
for(k in 1:7){
  data.graph1 <- NULL
  data.graph2 <- list(year16 = list(), year17 = list(), year18 = list(),
                      year19 = list(), year20 = list(), year21 = list()) 
  i=1
  for(year in 16:21){
    load(paste0("./PCA/pc.irlba/",pc[k],"pc/1000samples.",pc[k],"svd_year",year,".RData"))
    
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
                                Lower.genetic = Genetic.CI[,1],
                                Upper.genetic = Genetic.CI[,2],
                                True.genic.var = c(stage.year.Parents$genicG, stage.year.F1$genicG, stage.year.HDRW$genicG,
                                                   stage.year.PYT$genicG, stage.year.AYT$genicG, EYT1.true.Genic),
                                Estimated.genic.var = c(mean(Parents.Genic.var), mean(F1.Genic.var),
                                                        mean(HDRW.Genic.var), mean(PYT.Genic.var),
                                                        mean(AYT.Genic.var), mean(EYT1.Genic.var)),
                                Lower.genic = Genic.CI[,1],
                                Upper.genic = Genic.CI[,2])
    data.graph1 <- rbind(data.graph1, table_results)
    data.graph2[[i]] <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                                   PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                                   Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                                 PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                                   Stage = c(rep("Parents",1000), rep("F1",1000), rep("HDRW",1000),
                                             rep("PYT",1000), rep("AYT",1000), rep("EYT",1000)))
    i=i+1
  }
  if(pc[k] == 10){
    data_pc10.graph1 <- data.graph1
    data_pc10.graph2 <- data.graph2
  }else if(pc[k] == 50){
    data_pc50.graph1 <- data.graph1
    data_pc50.graph2 <- data.graph2
  }else if(pc[k] == 100){
    data_pc100.graph1 <- data.graph1
    data_pc100.graph2 <- data.graph2
  }else if(pc[k] == 500){
    data_pc500.graph1 <- data.graph1
    data_pc500.graph2 <- data.graph2
  }else if(pc[k] == 1000){
    data_pc1000.graph1 <- data.graph1
    data_pc1000.graph2 <- data.graph2
  }else  if(pc[k] == 2000){
    data_pc2000.graph1 <- data.graph1
    data_pc2000.graph2 <- data.graph2
  }else if(pc[k] == 3410){
    data_pc3410.graph1 <- data.graph1
    data_pc3410.graph2 <- data.graph2
  }else {
    print("ERROR")
  }
}


####################
# Graph: Bias
head(data_pc10.graph1)
str(data_pc10.graph2[1:2])

bias.data <- data.frame("Stage" = rep(data_pc10.graph1$Stage, 14))
bias.data$Year <- rep(rep(16:21, each=6), 14)
bias.data$Estimates <- rep(rep(c("Genetic", "Genic"), each=36), 7)
bias.data$pc.irlba <- rep(c("10", "50", "100", "500", "1000", "2000", "3410"), each=72)
bias.data$Value <- c((data_pc10.graph1$True.genetic.var - data_pc10.graph1$Estimated.genetic.var),
                     (data_pc10.graph1$True.genic.var - data_pc10.graph1$Estimated.genic.var),
                     (data_pc50.graph1$True.genetic.var - data_pc50.graph1$Estimated.genetic.var),
                     (data_pc50.graph1$True.genic.var - data_pc50.graph1$Estimated.genic.var),
                     (data_pc100.graph1$True.genetic.var - data_pc100.graph1$Estimated.genetic.var),
                     (data_pc100.graph1$True.genic.var - data_pc100.graph1$Estimated.genic.var),
                     (data_pc500.graph1$True.genetic.var - data_pc500.graph1$Estimated.genetic.var),
                     (data_pc500.graph1$True.genic.var - data_pc500.graph1$Estimated.genic.var),
                     (data_pc1000.graph1$True.genetic.var - data_pc1000.graph1$Estimated.genetic.var),
                     (data_pc1000.graph1$True.genic.var - data_pc1000.graph1$Estimated.genic.var),
                     (data_pc2000.graph1$True.genetic.var - data_pc2000.graph1$Estimated.genetic.var),
                     (data_pc2000.graph1$True.genic.var - data_pc2000.graph1$Estimated.genic.var),
                     (data_pc3410.graph1$True.genetic.var - data_pc3410.graph1$Estimated.genetic.var),
                     (data_pc3410.graph1$True.genic.var - data_pc3410.graph1$Estimated.genic.var))
bias.data$Stage <- factor(bias.data$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
bias.data$pc.irlba <- factor(bias.data$pc.irlba, levels = c("10", "50", "100", "500", "1000", "2000", "3410"))
head(bias.data)

optns <- theme_bw(base_size = 14, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_est <- c("#5FA83B", "#C11773")

pdf("./graph1_bias.pdf", width=20, height=13)
graph1 <- ggplot(bias.data, aes(y = Value, x = pc.irlba, group = Estimates)) +
  geom_line(aes(color=Estimates, linetype=Estimates)) +
  facet_grid(rows = vars(Stage), cols = vars(Year)) +
  scale_color_manual(values = color_est) +
  xlab("Number of Principal Components") + ylab("Bias") +
  geom_hline(aes(yintercept = 0), colour = "black", alpha=0.5, size=0.4, linetype="dashed")
graph1 + optns
dev.off()



