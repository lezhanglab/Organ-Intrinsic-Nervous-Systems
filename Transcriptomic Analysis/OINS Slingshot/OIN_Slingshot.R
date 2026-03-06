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

# get clusters using portal low-dim embedding
obj <- FindNeighbors(obj, reduction = "portal", dims = 1:20)
obj <- FindClusters(obj, resolution = .6)
DimPlot(obj, label = TRUE)
obj <- obj[, as.vector(obj@meta.data$seurat_clusters) != "12"]
obj <- FindNeighbors(obj, reduction = "portal", dims = 1:20)
# obj <- FindClusters(obj, resolution = 1.)
obj <- FindClusters(obj, resolution = .4)
DimPlot(obj, label = TRUE)


### Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(rownames(obj), row.names = rownames(obj))
colnames(gene_annotation) <- "gene_short_name"
# part two, cell information
cell_metadata <- as.data.frame(colnames(obj), row.names = colnames(obj))
colnames(cell_metadata) <- "barcode"
# part three, counts sparse matrix
New_matrix <- obj@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(obj), ]
expression_matrix <- New_matrix


### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


### Assign the cluster info
list_cluster <- obj@meta.data[["seurat_clusters"]]
names(list_cluster) <- obj@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate
# cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 100)
# cds_from_seurat <- reduce_dimension(cds_from_seurat)
# plot_cells(cds_from_seurat)
# cds_from_seurat@reducedDims@listData[["UMAP"]] <- obj@reductions[["umap"]]@cell.embeddings
# plot_cells(cds_from_seurat)
reducedDims(cds_from_seurat)$UMAP <- obj@reductions[["umap"]]@cell.embeddings
plot_cells(cds_from_seurat)


### Slingshot
DimPlot(obj, label = TRUE)
cds_from_seurat <- slingshot(data = cds_from_seurat, 
                             clusterLabels = obj$seurat_clusters,
                             end.clus = c("7","5","0","3"),
                             start.clus = c("1","4"))
cds_from_seurat


### Overall trajectories
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_1, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col="yellow", pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')


### invidiual trajectory
# 1
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_1, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
summary(cds_from_seurat$slingPseudotime_1)
# 2
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_2, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
summary(cds_from_seurat$slingPseudotime_2)
# 3
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_3, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
summary(cds_from_seurat$slingPseudotime_3)
# 4
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_4, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
summary(cds_from_seurat$slingPseudotime_4)
# 5
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(cds_from_seurat$slingPseudotime_5, breaks=100)]
plot(reducedDims(cds_from_seurat)$UMAP, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
summary(cds_from_seurat$slingPseudotime_5)

obj@meta.data$slingPseudotime_1 <- cds_from_seurat$slingPseudotime_1
obj@meta.data$slingPseudotime_2 <- cds_from_seurat$slingPseudotime_2
obj@meta.data$slingPseudotime_3 <- cds_from_seurat$slingPseudotime_3
obj@meta.data$slingPseudotime_4 <- cds_from_seurat$slingPseudotime_4
obj@meta.data$slingPseudotime_5 <- cds_from_seurat$slingPseudotime_5


### Min-max normalization for pseudo time / Quantile normalization
minmax_norm <- function(vec){
  vec_val <- vec[!is.na(vec)]
  vec_norm <- (vec_val - min(vec_val)) / (max(vec_val) - min(vec_val))
  vec[!is.na(vec)] <- vec_norm
  return(vec)
}
cds_from_seurat$slingPseudotime_1 <- minmax_norm(cds_from_seurat$slingPseudotime_1)
cds_from_seurat$slingPseudotime_2 <- minmax_norm(cds_from_seurat$slingPseudotime_2)
cds_from_seurat$slingPseudotime_3 <- minmax_norm(cds_from_seurat$slingPseudotime_3)
cds_from_seurat$slingPseudotime_4 <- minmax_norm(cds_from_seurat$slingPseudotime_4)
slingPseudotime(cds_from_seurat)
obj@meta.data$slingPseudotime_1_minimax_norm <- cds_from_seurat$slingPseudotime_1
obj@meta.data$slingPseudotime_2_minimax_norm <- cds_from_seurat$slingPseudotime_2
obj@meta.data$slingPseudotime_3_minimax_norm <- cds_from_seurat$slingPseudotime_3
obj@meta.data$slingPseudotime_4_minimax_norm <- cds_from_seurat$slingPseudotime_4
obj@meta.data$slingPseudotime_5_minimax_norm <- cds_from_seurat$slingPseudotime_5


### Save results
setwd("/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_trajectory_pseudo_time_analysis/results_slingshot_20240401")
saveRDS(obj, file = "obj_seurat.rds")
saveRDS(cds_from_seurat, file = "cds_slingshot.rds")
saveRDS(slingPseudotime(cds_from_seurat), file = "slingPseudotime.rds")


### Only keep slingPseudotime 1-4
# pst <- data.frame("Lineage1" = cds_from_seurat$slingPseudotime_1,
#                   "Lineage2" = cds_from_seurat$slingPseudotime_2,
#                   "Lineage3" = cds_from_seurat$slingPseudotime_3,
#                   "Lineage4" = cds_from_seurat$slingPseudotime_4)
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
pst <- pst[is.na(cds_from_seurat$slingPseudotime_5), ]
cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst
# plot
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(obj@meta.data$pseudo.time, breaks=100)]
plot(obj@reductions$umap@cell.embeddings, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
# 1
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
pst_tmp <- obj@meta.data$pseudo.time
pst_tmp[is.na(obj@meta.data$slingPseudotime_1_minimax_norm)] <- NA
plotcol <- colors[cut(pst_tmp, breaks=100)]
plot(obj@reductions$umap@cell.embeddings, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
# 2
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
pst_tmp <- obj@meta.data$pseudo.time
pst_tmp[is.na(obj@meta.data$slingPseudotime_2_minimax_norm)] <- NA
plotcol <- colors[cut(pst_tmp, breaks=100)]
plot(obj@reductions$umap@cell.embeddings, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
# 3
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
pst_tmp <- obj@meta.data$pseudo.time
pst_tmp[is.na(obj@meta.data$slingPseudotime_3_minimax_norm)] <- NA
plotcol <- colors[cut(pst_tmp, breaks=100)]
plot(obj@reductions$umap@cell.embeddings, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
# 4
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
pst_tmp <- obj@meta.data$pseudo.time
pst_tmp[is.na(obj@meta.data$slingPseudotime_4_minimax_norm)] <- NA
plotcol <- colors[cut(pst_tmp, breaks=100)]
plot(obj@reductions$umap@cell.embeddings, col=plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black')
# save obj
saveRDS(obj, file = "obj_pst.rds")
# save pst
pst_df <- obj@meta.data[, c("slingPseudotime_1", "slingPseudotime_2",
                            "slingPseudotime_3", "slingPseudotime_4",
                            "slingPseudotime_1_minimax_norm", "slingPseudotime_2_minimax_norm",
                            "slingPseudotime_3_minimax_norm", "slingPseudotime_4_minimax_norm",
                            "pseudo.time")]
write.csv(pst_df, file = "meta_pst.csv")

