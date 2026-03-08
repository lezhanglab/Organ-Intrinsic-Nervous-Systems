</> R
## R version 4.1.2 (2021-11-01)
## SeuratObject_4.1.3 sp_2.2-0    

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

####  Integrate ENS+Heart and ENS+Gut ####  
reference.list <- list(eh_2,
                       eg_1,
                       eg_2,
                       eh_1,
                       eg_4,
                       eg_3)
reference.list <- lapply(X = reference.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = reference.list)
Anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = features, dims = 1:50) 

coculture <- IntegrateData(anchorset = Anchors) 
coculture 
coculture <- ScaleData(object = coculture, verbose = FALSE)
coculture <- RunPCA(object = coculture, verbose = FALSE, npcs = 50)
DimPlot(object = coculture)

coculture <- RunUMAP(object = coculture, dims = 1:30) 
coculture <- FindNeighbors(object = coculture, dims = 1:30)
coculture <- FindClusters(object = coculture, resolution = 0.6)
head(coculture[[]])
DimPlot(object = coculture, reduction = 'umap', label = TRUE, shuffle = TRUE)
DimPlot(object = coculture, reduction = 'umap', label = TRUE, shuffle = TRUE, split.by = "orig.ident")

Idents(coculture) <- "seurat_clusters"
DimPlot(object = coculture, reduction = 'umap', label = TRUE, shuffle = TRUE)
g <- subset(coculture, orig.ident %in% c("ENS_Gut_042023","ENS_Gut_051923","ENS_Gut_081423",
                                   "ENS_Gut_090523"))
h <- subset(coculture, orig.ident %in% c("ENS_Heart_051923","ENS_Heart_090523"))
coculture$ID = coculture$orig.ident
coculture$ID[cells = WhichCells(g)]<- 'ENS_Gut' 
coculture$ID[cells = WhichCells(h)]<- 'ENS_Heart' 
DimPlot(object = coculture, reduction = 'umap', shuffle = TRUE, label = T)
DimPlot(object = coculture, reduction = 'umap', shuffle = TRUE, group.by = "ID", label = T)
DimPlot(object = coculture, reduction = 'umap', shuffle = TRUE, split.by = "ID", label = T)
saveRDS(coculture, "coculture.rds")

#####  Defining coculture states #####
## 1. Isolate state-defined genes using the 4-OINS integration dataset.
obj <- readRDS("4_OINS.rds")
Idents(obj) <- "seurat_clusters"
DimPlot(obj, label = T, shuffle = T)
pre <- subset(obj, seurat_clusters %in% c(1,4)) #precursor -> corresponding to Extended Data Fig. 14a, cluster 1,2
nb <- subset(obj, seurat_clusters %in% c(2,6)) #neuroblast -> corresponding to Extended Data Fig. 14a, cluster 3,4
n <- subset(obj, seurat_clusters %in% c(0,3,5,7)) #neuron -> corresponding to Extended Data Fig. 14a, cluster cA, cB, cC, cD
obj$state = obj$orig.ident
head(obj[[]])
obj$state[cells = WhichCells(pre)]<- 'Precursors'
obj$state[cells = WhichCells(nb)]<- 'Neuroblast'
obj$state[cells = WhichCells(n)]<- 'Neurons'
Idents(obj) <- "state"
DimPlot(obj, shuffle = T)

DefaultAssay(obj) <- "RNA"
obj_precursor_d <- FindMarkers(obj, ident.1 = 'Precursors', ident.2 = c("Neuroblast","Neurons"), only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox')
obj_precursor_d <- obj_precursor_d %>% arrange(desc(avg_log2FC))  
write.table(obj_precursor_d, "common_precursor_genes.xls",  sep = "\t")

obj_neuroblast_d <- FindMarkers(obj, ident.1 = 'Neuroblast', ident.2 = c("Precursors","Neurons"), only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox')
obj_neuroblast_d <- obj_neuroblast_d %>% arrange(desc(avg_log2FC))  
write.table(obj_neuroblast_d, "common_neuroblasts_genes.xls",  sep = "\t")

obj_neuron_d <- FindMarkers(obj, ident.1 = 'Neurons', ident.2 = c("Precursors","Neuroblast"), only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox')
obj_neuron_d <- obj_neuron_d %>% arrange(desc(avg_log2FC))  
write.table(obj_neuron_d, "common_neuron_genes.xls",  sep = "\t")

obj_precursor_d <- read.table("common_precursor_genes.xls",  sep = "\t")
obj_neuroblast_d <- read.table("common_neuroblasts_genes.xls",  sep = "\t")
obj_neuron_d <- read.table("common_neuron_genes.xls",  sep = "\t")

## Identify seurat_clusters of "coculture.rds" enriched for specific state genes, using high-resolution clustering
obj <- coculture
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(object = obj, resolution = 4) 
DimPlot(obj, shuffle = TRUE, label = TRUE) 
DefaultAssay(obj) <- "RNA"
Marker <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top <- Marker %>%
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 50) 

## Estimate the state-specific gene distributions in coculture clusters and define objects of "pre", "nb", "n", and "npc" 
neuron_cluster <- subset(top, top$gene %in% rownames(obj_neuron_d)) 
neuroblast_cluster <- subset(top, top$gene %in% rownames(obj_neuroblast_d)) 
precursor_cluster <- subset(top, top$gene %in% rownames(obj_precursor_d)) 


coculture$state = coculture$orig.ident
head(coculture[[]])
coculture$state[cells = WhichCells(pre)]<- 'Precursors'
coculture$state[cells = WhichCells(nb)]<- 'Neuroblast'
coculture$state[cells = WhichCells(n)]<- 'Neurons'
coculture$state[cells = WhichCells(npc)]<- 'NPC'

saveRDS(coculture, "coculture.rds")






