rm(list = ls())
library(ggplot2)
library(Seurat)
library(monocle3)
library(slingshot)
library(RColorBrewer)
library(scales)
# show_col(viridis_pal()(10))


### Trajectory 1-4
obj <- readRDS("obj_seurat.rds")
cds_from_seurat <- readRDS("cds_slingshot.rds")


### Only keep slingPseudotime 1-4
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
idx <- !(is.na(cds_from_seurat$slingPseudotime_1) & is.na(cds_from_seurat$slingPseudotime_2) & is.na(cds_from_seurat$slingPseudotime_3) & is.na(cds_from_seurat$slingPseudotime_4))
pst_bg <- pst[!idx, ]
pst <- pst[idx, ]

cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj_bg <- obj[, rownames(pst_bg)]
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst


pdf("slingshot_traj1to4_colorbyorgans.pdf", width = 5.5, height = 4.2)
plotcol <- as.vector(obj@meta.data$orig.ident)
plotcol[plotcol == "E14.5_Lung"] <- "#DE9BF4"
plotcol[plotcol == "E14.5_pancrease"] <- "#80E380"
plotcol[plotcol == "ICN.E14.3"] <- "#FF9380"
plotcol[plotcol == "ENS_E14"] <- "#8299FF"
plotcol_bg <- as.vector(obj_bg@meta.data$orig.ident)
plotcol_bg[plotcol_bg == "E14.5_Lung"] <- "#DE9BF4"
plotcol_bg[plotcol_bg == "E14.5_pancrease"] <- "#80E380"
plotcol_bg[plotcol_bg == "ICN.E14.3"] <- "#FF9380"
plotcol_bg[plotcol_bg == "ENS_E14"] <- "#8299FF"
plotcol <- c(plotcol_bg, plotcol)
umap_coor <- rbind(obj_bg@reductions$umap@cell.embeddings, obj@reductions$umap@cell.embeddings)
plot(umap_coor, col=plotcol, pch=16, asp = 1, cex=.3) + scale_fill_continuous(guide = "colourbar")
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black', linInd = 1:4)
dev.off()



### Trajectory 1
obj <- readRDS("obj_seurat.rds")
cds_from_seurat <- readRDS("cds_slingshot.rds")


### Only keep slingPseudotime 1-4
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
idx <- !(is.na(cds_from_seurat$slingPseudotime_1))
pst_bg <- pst[!idx, ]
pst <- pst[idx, ]

cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj_bg <- obj[, rownames(pst_bg)]
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst


pdf("slingshot_traj1_colorbyorgans.pdf", width = 5.5, height = 4.2)
plotcol <- as.vector(obj@meta.data$orig.ident)
plotcol[plotcol == "E14.5_Lung"] <- "#DE9BF4"
plotcol[plotcol == "E14.5_pancrease"] <- "#80E380"
plotcol[plotcol == "ICN.E14.3"] <- "#FF9380"
plotcol[plotcol == "ENS_E14"] <- "#8299FF"
plotcol_bg <- rep("#D3D3D3", dim(obj_bg)[2])
plotcol <- c(plotcol_bg, plotcol)
umap_coor <- rbind(obj_bg@reductions$umap@cell.embeddings, obj@reductions$umap@cell.embeddings)
plot(umap_coor, col=plotcol, pch=16, asp = 1, cex=.3) + scale_fill_continuous(guide = "colourbar")
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black', linInd = 1)
dev.off()



### Trajectory 2
obj <- readRDS("obj_seurat.rds")
cds_from_seurat <- readRDS("cds_slingshot.rds")


### Only keep slingPseudotime 1-4
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
idx <- !(is.na(cds_from_seurat$slingPseudotime_2))
pst_bg <- pst[!idx, ]
pst <- pst[idx, ]

cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj_bg <- obj[, rownames(pst_bg)]
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst


pdf("slingshot_traj2_colorbyorgans.pdf", width = 5.5, height = 4.2)
plotcol <- as.vector(obj@meta.data$orig.ident)
plotcol[plotcol == "E14.5_Lung"] <- "#DE9BF4"
plotcol[plotcol == "E14.5_pancrease"] <- "#80E380"
plotcol[plotcol == "ICN.E14.3"] <- "#FF9380"
plotcol[plotcol == "ENS_E14"] <- "#8299FF"
plotcol_bg <- rep("#D3D3D3", dim(obj_bg)[2])
plotcol <- c(plotcol_bg, plotcol)
umap_coor <- rbind(obj_bg@reductions$umap@cell.embeddings, obj@reductions$umap@cell.embeddings)
plot(umap_coor, col=plotcol, pch=16, asp = 1, cex=.3) + scale_fill_continuous(guide = "colourbar")
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black', linInd = 2)
dev.off()



### Trajectory 3
obj <- readRDS("obj_seurat.rds")
cds_from_seurat <- readRDS("cds_slingshot.rds")


### Only keep slingPseudotime 1-4
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
idx <- !(is.na(cds_from_seurat$slingPseudotime_3))
pst_bg <- pst[!idx, ]
pst <- pst[idx, ]

cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj_bg <- obj[, rownames(pst_bg)]
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst


pdf("slingshot_traj3_colorbyorgans.pdf", width = 5.5, height = 4.2)
plotcol <- as.vector(obj@meta.data$orig.ident)
plotcol[plotcol == "E14.5_Lung"] <- "#DE9BF4"
plotcol[plotcol == "E14.5_pancrease"] <- "#80E380"
plotcol[plotcol == "ICN.E14.3"] <- "#FF9380"
plotcol[plotcol == "ENS_E14"] <- "#8299FF"
plotcol_bg <- rep("#D3D3D3", dim(obj_bg)[2])
plotcol <- c(plotcol_bg, plotcol)
umap_coor <- rbind(obj_bg@reductions$umap@cell.embeddings, obj@reductions$umap@cell.embeddings)
plot(umap_coor, col=plotcol, pch=16, asp = 1, cex=.3) + scale_fill_continuous(guide = "colourbar")
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black', linInd = 3)
dev.off()



### Trajectory 4
obj <- readRDS("obj_seurat.rds")
cds_from_seurat <- readRDS("cds_slingshot.rds")


### Only keep slingPseudotime 1-4
pst <- slingPseudotime(cds_from_seurat)
pst <- pst[, c("Lineage1","Lineage2","Lineage3","Lineage4")]
rownames(pst) <- rownames(slingPseudotime(cds_from_seurat))
idx <- !(is.na(cds_from_seurat$slingPseudotime_4))
pst_bg <- pst[!idx, ]
pst <- pst[idx, ]

cal_mean <- function(vec){
  vec_val <- vec[!is.na(vec)]
  return(mean(vec_val))
}
pst <- apply(pst, 1, cal_mean)
obj_bg <- obj[, rownames(pst_bg)]
obj <- obj[, names(pst)]
obj@meta.data$pseudo.time <- pst


pdf("slingshot_traj4_colorbyorgans.pdf", width = 5.5, height = 4.2)
plotcol <- as.vector(obj@meta.data$orig.ident)
plotcol[plotcol == "E14.5_Lung"] <- "#DE9BF4"
plotcol[plotcol == "E14.5_pancrease"] <- "#80E380"
plotcol[plotcol == "ICN.E14.3"] <- "#FF9380"
plotcol[plotcol == "ENS_E14"] <- "#8299FF"
plotcol_bg <- rep("#D3D3D3", dim(obj_bg)[2])
plotcol <- c(plotcol_bg, plotcol)
umap_coor <- rbind(obj_bg@reductions$umap@cell.embeddings, obj@reductions$umap@cell.embeddings)
plot(umap_coor, col=plotcol, pch=16, asp = 1, cex=.3) + scale_fill_continuous(guide = "colourbar")
lines(SlingshotDataSet(cds_from_seurat), lwd=2, col='black', linInd = 4)
dev.off()


