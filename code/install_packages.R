#install the required packages

list.of.packages <- c("Seurat", "tidyr", "dplyr", "ggplot2", "cowplot", "gridExtra", 
                      "rjson", "harmony", "xlsx", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRepertoire")
