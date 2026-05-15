install.packages(c(
  "shiny",
  "plotly",
  "DT",
  "ggplot2",
  "dplyr",
  "tidyr",
  "networkD3",
  "RColorBrewer",
  "htmlwidgets"
))

if (!require("remotes")) install.packages("remotes")
remotes::install_github("scverse/anndataR")

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "rhdf5"))