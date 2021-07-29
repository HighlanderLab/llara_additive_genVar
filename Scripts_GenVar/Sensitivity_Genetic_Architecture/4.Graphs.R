
#######################################################################
#          Genetic and genic variances obtained from julia            #
#           Letícia Lara, Thiago Oliveira, Gregor Gorjanc             #
#                            June, 2021                               #
#######################################################################

rm(list = ls())
library(ggplot2)             # Version 3.3.3
library(coda)                # Version 0.19.4
library(tidyr)               # Version 1.1.2


#### Loading and manipulating the data ####
QTLs = c(10, 100, 1000)
nSNPs = 500
scenGxY = c("noGxY", "GxY")
dir = "alphabayes"
startYear = 16
Yrs = 6
nSamples = 1000


#### Running all repetitions ####
for(r in 1:10){
  for(Scen in scenGxY){
    for(nQTL in QTLs){
      load(paste0("./AlphaSimR/sim", Sim, "/Burnin_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
      load(paste0("./AlphaSimR/sim", Sim, "/all_stages_and_years_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
      
      
      ####################
      ### Sumarizing results - Table
      data.graph1 <- NULL
      data.graph2 <- list(year16 = list(), year17 = list(), year18 = list(),
                          year19 = list(), year20 = list(), year21 = list())
      i=1
      
      for(year in sizeTrain:21){
        load(paste0("./", dir, "/estimates/sim", Sim, "/Estimates_", dir, "_year", year, "_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
        
        Genetic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genetic.var), prob = 0.95), 
                            HPDinterval(as.mcmc(F1.Genetic.var), prob = 0.95), 
                            HPDinterval(as.mcmc(HDRW.Genetic.var), prob = 0.95), 
                            HPDinterval(as.mcmc(PYT.Genetic.var), prob = 0.95), 
                            HPDinterval(as.mcmc(AYT.Genetic.var), prob = 0.95), 
                            HPDinterval(as.mcmc(EYT1.Genetic.var), prob = 0.95))
        Genic.CI <- rbind(HPDinterval(as.mcmc(Parents.Genic.var), prob = 0.95), 
                          HPDinterval(as.mcmc(F1.Genic.var), prob = 0.95), 
                          HPDinterval(as.mcmc(HDRW.Genic.var), prob = 0.95), 
                          HPDinterval(as.mcmc(PYT.Genic.var), prob = 0.95), 
                          HPDinterval(as.mcmc(AYT.Genic.var), prob = 0.95), 
                          HPDinterval(as.mcmc(EYT1.Genic.var), prob = 0.95))
        
        table_results <- data.frame(Stage = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"),
                                    True.genetic.var = c(Parents.true.Genetic, F1.true.Genetic, HDRW.true.Genetic,
                                                         PYT.true.Genetic, AYT.true.Genetic, EYT1.true.Genetic),
                                    Estimated.genetic.var = c(mean(Parents.Genetic.var, na.rm = TRUE), mean(F1.Genetic.var, na.rm = TRUE),
                                                              mean(HDRW.Genetic.var, na.rm = TRUE), mean(PYT.Genetic.var, na.rm = TRUE),
                                                              mean(AYT.Genetic.var, na.rm = TRUE), mean(EYT1.Genetic.var, na.rm = TRUE)),
                                    Lower.genetic = Genetic.CI[,1],
                                    Upper.genetic = Genetic.CI[,2],
                                    True.genic.var = c(Parents.true.Genic, F1.true.Genic, HDRW.true.Genic,
                                                       PYT.true.Genic, AYT.true.Genic, EYT1.true.Genic),
                                    Estimated.genic.var = c(mean(Parents.Genic.var, na.rm = TRUE), mean(F1.Genic.var, na.rm = TRUE),
                                                            mean(HDRW.Genic.var, na.rm = TRUE), mean(PYT.Genic.var, na.rm = TRUE),
                                                            mean(AYT.Genic.var, na.rm = TRUE), mean(EYT1.Genic.var, na.rm = TRUE)),
                                    Lower.genic = Genic.CI[,1],
                                    Upper.genic = Genic.CI[,2],
                                    True.CovG_HW = c(Parents.CovG_HW.true, F1.CovG_HW.true, HDRW.CovG_HW.true,
                                                     PYT.CovG_HW.true, AYT.CovG_HW.true, EYT1.CovG_HW.true))
        
        data.graph1 <- rbind(data.graph1, table_results)
        data.graph2[[i]] <- data.frame(Genetic.var = c(Parents.Genetic.var, F1.Genetic.var, HDRW.Genetic.var,
                                                       PYT.Genetic.var, AYT.Genetic.var, EYT1.Genetic.var),
                                       Genic.var = c(Parents.Genic.var, F1.Genic.var, HDRW.Genic.var,
                                                     PYT.Genic.var, AYT.Genic.var, EYT1.Genic.var),
                                       Stage = c(rep("Parents", nSamples), rep("F1", nSamples), rep("HDRW", nSamples),
                                                 rep("PYT", nSamples), rep("AYT", nSamples), rep("EYT", nSamples)))
        data.graph2[[i]]$Stage <- factor(data.graph2[[i]]$Stage, levels = c("EYT", "AYT", "PYT", "HDRW", "F1", "Parents"))
        i = i + 1
      }
      save(data.graph1, data.graph2, file = paste0("./", dir, "/estimates/sim", Sim, "/Data.graph_", dir, "_", Scen, "_nQTLs", nQTL, "_Sim", Sim, ".RData"))
    }
  }
}



##### Graph: Estimated vs True variances - Scatter Plot  #####

rm(list = ls())
dir = "alphabayes"


# Function: Loading data and preparing graphs
scatter.graphs <- function(GxY = FALSE, nQTL, var = "genetic", method = dir){
  if(GxY == TRUE){
    intGxY = "GxY"
  } else{
    intGxY = "noGxY"
  }
  shape_years <- c(16, 17, 18, 7, 3, 8)
  
  graph3.genetic_final = graph3.genic_final = NULL
  for(nSim in 1:10){
    load(paste0("./", method, "/sim", nSim, "/Data.graph_", method, "_", intGxY, "_nQTLs", nQTL, "_Sim", nSim, ".RData"))
    data.graph3 <- as.data.frame(data.graph1)
    data.graph3 = mutate(data.graph3, Genic.CovG = True.genic.var + True.CovG_HW)
    data.graph3$Year <- rep(16:21, each = 6)
    data.graph3$Stage <- factor(data.graph3$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
    
    graph3.genetic <- data.graph3[, c("Stage", "Year", "True.genetic.var", "Estimated.genetic.var")]
    graph3.genetic$Stage <- factor(graph3.genetic$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
    graph3.genetic$Year <- as.factor(graph3.genetic$Year)
    graph3.genetic$Sim = as.factor(nSim)
    graph3.genic <- data.graph3[, c("Stage", "Year", "Genic.CovG", "Estimated.genic.var")]
    graph3.genic$Stage <- factor(graph3.genic$Stage, levels = c("Parents", "F1", "HDRW", "PYT", "AYT", "EYT"))
    graph3.genic$Year <- as.factor(graph3.genic$Year)
    graph3.genic$Sim = as.factor(nSim)
    
    graph3.genetic_final = rbind(graph3.genetic_final, graph3.genetic)
    graph3.genic_final = rbind(graph3.genic_final, graph3.genic)
  }
  
  optns <- theme_bw(base_size = 10, base_family = "sans") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  update_geom_defaults("point", list(size = 1))
  theme_set(theme_grey(base_size = 2))
  
  color_rgb <- c("#666666", "#000000", "#C11773", "#5FA83B", "#DB4F0F", "#3798D9")
  color_ci <- c("#E5E5E5", "#666666", "#DE99BF", "#B2E39C", "#DFB19C", "#9DC4E0")
  
  if(var == "genetic"){
    graph3 <- ggplot(graph3.genetic_final, aes(True.genetic.var, Estimated.genetic.var)) + 
      geom_point(aes(color = Stage, shape = Year), size = 2) + coord_cartesian(ylim = c(0, 0.85), xlim = c(0, 0.25)) +
      geom_abline(aes(intercept = 0, slope = 1), colour = "red", alpha = 0.5, size = 0.4, linetype = "dashed") +
      labs(x = "True genetic variance", y = "Estimated genetic variance") +
      scale_color_manual(values = color_rgb) + 
      scale_shape_manual(values = shape_years) + optns + guides(fill = FALSE)
    result = graph3
  } 
  else if(var == "genic"){
    graph3 <- ggplot(graph3.genic_final, aes(Genic.CovG, Estimated.genic.var)) + 
      geom_point(aes(color = Stage, shape = Year), size = 2) + coord_cartesian(ylim = c(0, 0.85), xlim = c(0, 0.25)) +
      geom_abline(aes(intercept = 0, slope = 1), colour = "red", alpha = 0.5, size = 0.4, linetype = "dashed") +
      labs(x = "True genic variance", y = "Estimated genic variance") +
      scale_color_manual(values = color_rgb) + 
      scale_shape_manual(values = shape_years) + optns + guides(fill = FALSE)
    result = graph3
  }  
  return(result)
}

# Obtaining graphs - Genetic variance
g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)}
graph_legend = scatter.graphs(GxY = FALSE, nQTL = 10, var = "genetic")
mylegend = g_legend(graph_legend + theme(legend.position="bottom"))

pdf(paste0("./scatter_Genetic_", dir, ".pdf"), width = 12, height = 18)
graph_noGxY_nQTLs10 = scatter.graphs(GxY = FALSE, nQTL = 10, var = "genetic")
graph_noGxY_nQTLs100 = scatter.graphs(GxY = FALSE, nQTL = 100, var = "genetic")
graph_noGxY_nQTLs1000 = scatter.graphs(GxY = FALSE, nQTL = 1000, var = "genetic")
graph_GxY_nQTLs10 = scatter.graphs(GxY = TRUE, nQTL = 10, var = "genetic")
graph_GxY_nQTLs100 = scatter.graphs(GxY = TRUE, nQTL = 100, var = "genetic")
graph_GxY_nQTLs1000 = scatter.graphs(GxY = TRUE, nQTL = 1000, var = "genetic")

grid.arrange(arrangeGrob(graph_noGxY_nQTLs10 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 10 QTLs")), 
                         graph_GxY_nQTLs10 + theme(legend.position = "none") + ggtitle(paste0("GxY and 10 QTLs")), 
                         graph_noGxY_nQTLs100 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 100 QTLs")), 
                         graph_GxY_nQTLs100 + theme(legend.position = "none") + ggtitle(paste0("GxY and 100 QTLs")), 
                         graph_noGxY_nQTLs1000 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 1000 QTLs")), 
                         graph_GxY_nQTLs1000 + theme(legend.position = "none") + ggtitle(paste0("GxY and 1000 QTLs")),
                         ncol = 2, nrow = 3),
             arrangeGrob(mylegend), heights = c(18, 0.8)) 
dev.off()

# Obtaining graphs - Genic variance
g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)}
graph_legend = scatter.graphs(GxY = FALSE, nQTL = 10, var = "genic")
mylegend = g_legend(graph_legend + theme(legend.position="bottom"))

pdf(paste0("./scatter_Genic_", dir, ".pdf"), width = 12, height = 18)
graph_noGxY_nQTLs10 = scatter.graphs(GxY = FALSE, nQTL = 10, var = "genic")
graph_noGxY_nQTLs100 = scatter.graphs(GxY = FALSE, nQTL = 100, var = "genic")
graph_noGxY_nQTLs1000 = scatter.graphs(GxY = FALSE, nQTL = 1000, var = "genic")
graph_GxY_nQTLs10 = scatter.graphs(GxY = TRUE, nQTL = 10, var = "genic")
graph_GxY_nQTLs100 = scatter.graphs(GxY = TRUE, nQTL = 100, var = "genic")
graph_GxY_nQTLs1000 = scatter.graphs(GxY = TRUE, nQTL = 1000, var = "genic")

grid.arrange(arrangeGrob(graph_noGxY_nQTLs10 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 10 QTLs")), 
                         graph_GxY_nQTLs10 + theme(legend.position = "none") + ggtitle(paste0("GxY and 10 QTLs")), 
                         graph_noGxY_nQTLs100 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 100 QTLs")), 
                         graph_GxY_nQTLs100 + theme(legend.position = "none") + ggtitle(paste0("GxY and 100 QTLs")), 
                         graph_noGxY_nQTLs1000 + theme(legend.position = "none") + ggtitle(paste0("noGxY and 1000 QTLs")), 
                         graph_GxY_nQTLs1000 + theme(legend.position = "none") + ggtitle(paste0("GxY and 1000 QTLs")),
                         ncol = 2, nrow = 3),
             arrangeGrob(mylegend), heights = c(18, 0.8)) 
dev.off()


