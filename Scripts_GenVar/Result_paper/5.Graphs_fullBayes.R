
#######################################################################
#               Creating Graphs from Gen Var analysis                 #
#                 Simulated wheat breeding program                    #
#                   Leticia Lara - February, 2020                     #
#######################################################################

library(ggplot2)
library(coda)
library(openxlsx)
library(tidyr)
library(ggridges)
library(ggpubr)
library(AlphaSimR)           # Version 0.11.0

rm(list=ls())

load("../Burnin.RData")
load("../all_stages_and_year.RData")

####################
### Sumarizing results - Table
data.graph1 <- NULL
data.graph2 <- list(year16 = list(), year17 = list(), year18 = list(),
                    year19 = list(), year20 = list(), year21 = list()) 
i=1

for(year in 16:21){
  load(paste0("./1000samples.fullBayes_year",year,".RData"))

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
  if(year == 16){
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = paste0("year_",year)) 
    writeData(wb, sheet = paste0("year_",year), x = table_results)
  }else{
    addWorksheet(wb, sheetName = paste0("year_",year)) 
    writeData(wb, sheet = paste0("year_",year), x = table_results)
  }
  
  data.graph1 <- rbind(data.graph1, table_results)
  data.graph2[[i]] <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                                 PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                                 Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                               PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                                 Stage = c(rep("Parents",1000), rep("F1",1000), rep("HDRW",1000),
                                           rep("PYT",1000), rep("AYT",1000), rep("EYT",1000)))
  data.graph2[[i]]$Stage <- factor(data.graph2[[i]]$Stage, levels = c("EYT", "AYT", "PYT", "HDRW", "F1", "Parents"))
  i=i+1
}
#saveWorkbook(wb, "table_results.xlsx")
#save(data.graph1, data.graph2, file="data.graphs.RData")
load("data.graphs.RData")


####################
# Graph 1: Genic and Genetic variances x Year for each stage
data.graph1 <- as.data.frame(data.graph1)
data.graph1$Year <- c(rep(16,6), rep(17,6), rep(18,6), rep(19,6), rep(20,6), rep(21,6))
data.graph1$Stage <- factor(data.graph1$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
head(data.graph1)

optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")

pdf("./graph1_genetic_1000SIM.pdf", width=11, height=10)
graph1.genetic <- gather(data.graph1, Genetic.var, Value, True.genetic.var:Estimated.genetic.var, factor_key = TRUE)
graph1.genetic$Stage <- factor(graph1.genetic$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
levels(graph1.genetic$Genetic.var) <- sub("True.genetic.var", "True", levels(graph1.genetic$Genetic.var))
levels(graph1.genetic$Genetic.var) <- sub("Estimated.genetic.var", "Estimated", levels(graph1.genetic$Genetic.var))
head(graph1.genetic)

graph1.1 <- ggplot(graph1.genetic, aes(y = Value, x = Year, group = Genetic.var)) +
  geom_ribbon(aes(ymin = Lower.genetic, ymax = Upper.genetic, fill = Stage), alpha=.2) +
  geom_line(aes(color = Stage, linetype = Genetic.var)) + 
  xlab("Year") + ylab("Genetic variance") +
  guides(color = FALSE) + facet_wrap(~ Stage) + labs(linetype = "Genetic variance")+
  scale_color_manual(values = color_rgb) +
  scale_fill_manual(values = color_ci) +
  coord_cartesian(ylim = c(0, 0.15)) + optns + guides(fill = FALSE)
graph1.1
dev.off()

pdf("./graph1_genic_1000SIM.pdf", width = 11, height = 10)
graph1.genic <- gather(data.graph1, Genic.var, Value, True.genic.var:Estimated.genic.var, factor_key = TRUE)
graph1.genic$Stage <- factor(graph1.genic$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
levels(graph1.genic$Genic.var) <- sub("True.genic.var", "True", levels(graph1.genic$Genic.var))
levels(graph1.genic$Genic.var) <- sub("Estimated.genic.var", "Estimated", levels(graph1.genic$Genic.var))
head(graph1.genic)

graph1.2 <- ggplot(graph1.genic, aes(y = Value, x = Year, group = Genic.var)) +
  geom_ribbon(aes(ymin = Lower.genic, ymax = Upper.genic, fill = Stage), alpha=.2) +
  geom_line(aes(color = Stage, linetype = Genic.var)) + 
  xlab("Year") + ylab("Genic variance") +
  guides(color = FALSE) + facet_wrap(~ Stage) + labs(linetype = "Genic variance")+
  scale_color_manual(values = color_rgb) +
  scale_fill_manual(values = color_ci) +
  coord_cartesian(ylim = c(0, 0.15)) + optns + guides(fill = FALSE)
graph1.2
dev.off()


####################
# Graph 2: Genetic variance distribution
head(data.graph2$year16)

true.genetic <- as.data.frame(c(data.graph1[(data.graph1$Stage == "Parents" & data.graph1$Year == 16), 2],
                                data.graph1[(data.graph1$Stage == "F1" & data.graph1$Year == 17), 2],
                                data.graph1[(data.graph1$Stage == "HDRW" & data.graph1$Year == 18), 2],
                                data.graph1[(data.graph1$Stage == "PYT" & data.graph1$Year == 19), 2],
                                data.graph1[(data.graph1$Stage == "AYT" & data.graph1$Year == 20), 2],
                                data.graph1[(data.graph1$Stage == "EYT" & data.graph1$Year == 21), 2]))
colnames(true.genetic) <- c("true.genetic"); true.genetic

true.genic <- as.data.frame(c(data.graph1[(data.graph1$Stage == "Parents" & data.graph1$Year == 16), 6],
                              data.graph1[(data.graph1$Stage == "F1" & data.graph1$Year == 17), 6],
                              data.graph1[(data.graph1$Stage == "HDRW" & data.graph1$Year == 18), 6],
                              data.graph1[(data.graph1$Stage == "PYT" & data.graph1$Year == 19), 6],
                              data.graph1[(data.graph1$Stage == "AYT" & data.graph1$Year == 20), 6],
                              data.graph1[(data.graph1$Stage == "EYT" & data.graph1$Year == 21), 6]))
colnames(true.genic) <- c("true.genic"); true.genic

new.data.graph2 <- as.data.frame(rbind(data.graph2$year16[data.graph2$year16$Stage == "Parents",],
                                       data.graph2$year17[data.graph2$year17$Stage == "F1",],
                                       data.graph2$year18[data.graph2$year18$Stage == "HDRW",],
                                       data.graph2$year19[data.graph2$year19$Stage == "PYT",],
                                       data.graph2$year20[data.graph2$year20$Stage == "AYT",],
                                       data.graph2$year21[data.graph2$year21$Stage == "EYT",]))
head(new.data.graph2)

optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#3798D9", "#DB4F0F", "#5FA83B", "#C11773", "#000000", "#666666")
color_ci <- c("#9DC4E0", "#DFB19C", "#B2E39C", "#DE99BF", "#666666", "#E5E5E5")

pdf("./graph2_genetic_distribution_1000SIM.pdf", width = 11, height = 10)
graph2.1 <- ggplot(new.data.graph2, aes(x = Genetic.var, y = Stage, fill = Stage)) +
  geom_density_ridges(alpha = 0.5, color = "white",
                      scale = 2, rel_min_height = 0.001) +
  labs(x = "Genetic variance", y = "Stage") + xlim(0, 0.125) +
  guides(fill = FALSE) + theme_ridges() +
  geom_segment(aes(x=true.genetic[1,], xend=true.genetic[1,], y = 5.9, yend = 6.2), size=1.2) +
  geom_segment(aes(x=true.genetic[2,], xend=true.genetic[2,], y = 4.9, yend = 5.2), size=1.2) +
  geom_segment(aes(x=true.genetic[3,], xend=true.genetic[3,], y = 3.9, yend = 4.2), size=1.2) +
  geom_segment(aes(x=true.genetic[4,], xend=true.genetic[4,], y = 2.9, yend = 3.2), size=1.2) +
  geom_segment(aes(x=true.genetic[5,], xend=true.genetic[5,], y = 1.9, yend = 2.2), size=1.2) +
  geom_segment(aes(x=true.genetic[6,], xend=true.genetic[6,], y = 0.9, yend = 1.2), size=1.2) +
  scale_size_identity() +
  scale_fill_manual(values = color_ci)
final.graph2.1 <- graph2.1 + theme_classic() + optns
suppressMessages(print(final.graph2.1))
dev.off()

pdf("./graph2_genic_distribution_1000SIM.pdf", width = 11, height = 10)
graph2.2 <- ggplot(new.data.graph2, aes(x = Genic.var, y = Stage, fill = Stage)) +
  geom_density_ridges(alpha = 0.5, color = "white",
                      scale = 2, rel_min_height = 0.001) +
  labs(x = "Genic variance", y = "Stage") + xlim(0, 0.125) +
  guides(fill = FALSE) + theme_ridges() +
  geom_segment(aes(x=true.genic[1,], xend=true.genic[1,], y = 5.9, yend = 6.2), size=1.2) +
  geom_segment(aes(x=true.genic[2,], xend=true.genic[2,], y = 4.9, yend = 5.2), size=1.2) +
  geom_segment(aes(x=true.genic[3,], xend=true.genic[3,], y = 3.9, yend = 4.2), size=1.2) +
  geom_segment(aes(x=true.genic[4,], xend=true.genic[4,], y = 2.9, yend = 3.2), size=1.2) +
  geom_segment(aes(x=true.genic[5,], xend=true.genic[5,], y = 1.9, yend = 2.2), size=1.2) +
  geom_segment(aes(x=true.genic[6,], xend=true.genic[6,], y = 0.9, yend = 1.2), size=1.2) +
  scale_size_identity() +
  scale_fill_manual(values = color_ci)
final.graph2.2 <- graph2.2 + theme_classic() + optns
suppressMessages(print(final.graph2.2))
dev.off()


####################
# Graph 3: True and estimated values for genetic and genic variance
data.graph3 <- data.graph1[,c(1, 2, 6, 3, 7, 10)]
head(data.graph3)

optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")

pdf("./graph3_true_1000SIM.pdf", width = 11, height = 10)
graph3.true <- gather(data.graph3, True.value, Value, True.genetic.var:True.genic.var, factor_key = TRUE)
graph3.true$Stage <- factor(graph3.true$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
levels(graph3.true$True.value) <- sub("True.genetic.var", "Genetic", levels(graph3.true$True.value))
levels(graph3.true$True.value) <- sub("True.genic.var", "Genic", levels(graph3.true$True.value))
graph3.true$True.value <- factor(graph3.true$True.value, levels = c("Genic", "Genetic"))
head(graph3.true)

graph3.1 <- ggplot(graph3.true, aes(y = Value, x = Year, group = True.value)) +
  geom_line(aes(color = Stage, linetype = True.value)) + 
  xlab("Year") + ylab("True value") +
  guides(color = FALSE) + facet_wrap(~ Stage) + labs(linetype = "True value")+
  scale_color_manual(values = color_rgb) +
  scale_fill_manual(values = color_ci) +
  coord_cartesian(ylim = c(0, 0.1))
graph3.1 + optns + guides(fill = FALSE)
dev.off()

pdf("./graph3_estimated_1000SIM.pdf", width = 11, height = 10)
graph3.estimated <- gather(data.graph3, Estimated.value, Value, Estimated.genetic.var:Estimated.genic.var, factor_key = TRUE)
graph3.estimated$Stage <- factor(graph3.estimated$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
levels(graph3.estimated$Estimated.value) <- sub("Estimated.genetic.var", "Genetic", levels(graph3.estimated$Estimated.value))
levels(graph3.estimated$Estimated.value) <- sub("Estimated.genic.var", "Genic", levels(graph3.estimated$Estimated.value))
graph3.estimated$Estimated.value <- factor(graph3.estimated$Estimated.value, levels = c("Genic", "Genetic"))
head(graph3.estimated)

graph3.1 <- ggplot(graph3.estimated, aes(y = Value, x = Year, group = Estimated.value)) +
  geom_line(aes(color = Stage, linetype = Estimated.value)) + 
  xlab("Year") + ylab("Estimated value") +
  guides(color = FALSE) + facet_wrap(~ Stage) + labs(linetype = "Estimated value")+
  scale_color_manual(values = color_rgb) +
  scale_fill_manual(values = color_ci) +
  coord_cartesian(ylim = c(0, 0.1))
graph3.1 + optns + guides(fill = FALSE)
dev.off()


####################
# Graph 4: True and estimated values for genetic and genic variance in the same graph
data.graph4 <- gather(data.graph3, Variance, Value, True.genetic.var:Estimated.genic.var, factor_key = TRUE)
data.graph4$Stage <- factor(data.graph4$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
data.graph4$Type <- c(rep("True",72), rep("Estimated",72))
levels(data.graph4$Variance) <- sub("Estimated.genetic.var", "Genetic", levels(data.graph4$Variance))
levels(data.graph4$Variance) <- sub("Estimated.genic.var", "Genic", levels(data.graph4$Variance))
levels(data.graph4$Variance) <- sub("True.genetic.var", "Genetic", levels(data.graph4$Variance))
levels(data.graph4$Variance) <- sub("True.genic.var", "Genic", levels(data.graph4$Variance))
data.graph4$Lower <- rep(c(data.graph1[,4], data.graph1[,8]), 2)
data.graph4$Upper <- rep(c(data.graph1[,5], data.graph1[,9]), 2)
data.graph4$Type <- factor(data.graph4$Type, levels = c("True", "Estimated"))
head(data.graph4)

optns <- theme_bw(base_size = 20, base_family = "sans") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=6))

color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")

pdf("./graph4.pdf", width = 11, height = 10)
graph4 <- ggplot(data.graph4, aes(y = Value, x = Year, linetype = Type, shape = Variance)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Stage), alpha=.2) +
  geom_line(aes(color = Stage)) + 
  geom_point(aes(color = Stage), size = 1.5) +
  xlab("Year") + ylab("Variance Value") +
  guides(color = FALSE) + facet_wrap(~ Stage) + labs(linetype = " ", shape = " ")+
  scale_color_manual(values = color_rgb) +
  scale_fill_manual(values = color_ci) +
  scale_shape_manual("", values=c(5,17))
graph4 + optns + guides(fill = FALSE)
dev.off()


##############################################################################################################
##############################################################################################################

