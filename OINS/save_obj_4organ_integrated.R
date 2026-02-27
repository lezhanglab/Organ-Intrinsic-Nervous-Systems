rm(list = ls())
library(ggplot2)
library(Seurat)
library(monocle3)
library(slingshot)
library(RColorBrewer)


### Reading in Seurat object
setwd("/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_trajectory_integration_portal_harmony/results_portal")
obj <- readRDS("Mar25_portal_4organOIN_withoutNPC_v1.rds")


### Run clustering
# store low-dim embedding from Portal
lowdim <- read.csv("Mar25_portal_4organOIN_withoutNPC_v1_lowdim.csv", row.names = 1)
get_barcode <- function(s){
  return(paste(strsplit(s, split = "-")[[1]][1], strsplit(s, split = "-")[[1]][2], sep = "-"))
}
rownames(lowdim) <- sapply(rownames(lowdim), get_barcode)
lowdim <- lowdim[colnames(obj), ]
lowdim <- as.matrix(lowdim)
reduction.save <- "portal"
reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
obj[[reduction.save]] <- Seurat::CreateDimReducObject(
  embeddings = lowdim,
  stdev = as.numeric(apply(lowdim, 2, stats::sd)),
  assay = "RNA",
  key = reduction.key
)


### load pseudotime info
meta_all <- read.csv("/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/ref_fig6_embryopaper/meta_all.csv", row.names = 1)
obj <- obj[, rownames(meta_all)]
DimPlot(obj, group.by = "orig.ident")
saveRDS(obj, file = "/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_trajectory_integration_portal_harmony/obj_4organ_integrated.rds")

