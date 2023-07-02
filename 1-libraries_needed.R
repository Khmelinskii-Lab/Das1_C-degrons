##### for project : Analysis of degron at C terminal
##### sub text:all libraries needed
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

# libraries needed 
library(dplyr)
library(ggplot2)
library(dummies)
library(ggpubr)
library(cowplot)
library(factoextra)
library(stringr)
library(stringi)
library(protr)
library(ggrepel)
library(venn)
library(iml)
library(readxl)
library(keras)
library(tensorflow)
install_keras()
library(Peptides)
library(reshape2)
library(limma)
library("ggpubr")
library(rstatix)



# create all directory

dir.create(paste0(path_plot,"miscelleneous/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

dir.create(paste0(path_plot,"Figure1/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

dir.create(paste0(path_plot,"Figure2/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)
dir.create(paste0(path_plot,"Figure3/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)
dir.create(
  paste0(path_plot,"Figure4/"),
  showWarnings = TRUE,
  recursive = FALSE,
  mode = "0777"
)

dir.create(
  paste0(path_plot,"Figure5/"),
  showWarnings = TRUE,
  recursive = FALSE,
  mode = "0777"
)

dir.create(paste0(path_plot,"Figure6/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)



dir.create(paste0(path_plot,"FigureS5/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)
dir.create(paste0(path_plot,"FigureS6/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

dir.create(paste0(path_plot,"FigureS1/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)
dir.create(
  paste0(path_plot,"FigureS2/"),
  showWarnings = TRUE,
  recursive = FALSE,
  mode = "0777"
)
dir.create(
  paste0(path_plot,"FigureS2/cluster_instable_12x/"),
  showWarnings = TRUE,
  recursive = FALSE,
  mode = "0777"
)
