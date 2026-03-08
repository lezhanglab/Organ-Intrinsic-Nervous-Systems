</> R
## R version 4.1.2 (2021-11-01)
## SeuratObject_4.1.3 sp_2.2-0    

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

#ENS+Heart 051923 and 090523
ENS_Heart_051923 <- Read10X(data.dir = "ENS_Heart_051923")
ENS_Heart_051923 <- CreateSeuratObject(counts = ENS_Heart_051923, project = "ENS_Heart_051923")
ENS_Heart_051923
#5146 cells
ENS_Heart_090523 <- Read10X(data.dir = "ENS_Heart_090523")
ENS_Heart_090523 <- CreateSeuratObject(counts = ENS_Heart_090523, project = "ENS_Heart_090523")
ENS_Heart_090523
#8994 cells

## Quanlity control 
ENS_Heart_051923
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Heart_051923), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Heart_051923, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Heart_051923, slot = 'counts'))
ENS_Heart_051923[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Heart_051923, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Heart_051923.QC <- subset(x = ENS_Heart_051923, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Heart_051923.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Heart_051923.QC #5063

ENS_Heart_090523 #8994
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Heart_090523), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Heart_090523, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Heart_090523, slot = 'counts'))
ENS_Heart_090523[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Heart_090523, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Heart_090523.QC <- subset(x = ENS_Heart_090523, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Heart_090523.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Heart_090523.QC #8521

## rename datasets, add sample id
ENS_Heart_051923.QC <- RenameCells(ENS_Heart_051923.QC, add.cell.id = "ENS_Heart_051923")
Cells(ENS_Heart_051923.QC)
ENS_Heart_051923.QC

ENS_Heart_090523.QC <- RenameCells(ENS_Heart_090523.QC, add.cell.id = "ENS_Heart_090523")
Cells(ENS_Heart_090523.QC)
ENS_Heart_090523.QC

#Unsupervised clustering
ENS_Heart_051923.QC.1 <- NormalizeData(object = ENS_Heart_051923.QC, verbose = FALSE)
ENS_Heart_051923.QC.1 <- FindVariableFeatures(object = ENS_Heart_051923.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Heart_051923.QC.1))
ENS_Heart_051923.QC.1 <- ScaleData(object = ENS_Heart_051923.QC.1, verbose = FALSE)
ENS_Heart_051923.QC.1 <- RunPCA(object = ENS_Heart_051923.QC.1, verbose = FALSE)
ENS_Heart_051923.QC.1 <- RunUMAP(object = ENS_Heart_051923.QC.1, dims = 1:50)
ENS_Heart_051923.QC.1 <- FindNeighbors(object = ENS_Heart_051923.QC.1, dims = 1:50)
ENS_Heart_051923.QC.1 <- FindClusters(object = ENS_Heart_051923.QC.1, resolution = 0.6)
DefaultAssay(object = ENS_Heart_051923.QC.1) <- "RNA"
VlnPlot(object = ENS_Heart_051923.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 2,5,6 and 15 are Phox2b+
VlnPlot(object = ENS_Heart_051923.QC.1, features = c("Sox10"), pt.size = 0.1) #None
VlnPlot(object = ENS_Heart_051923.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 2,5,6 and 15 are Snap25+
VlnPlot(object = ENS_Heart_051923.QC.1, features = c("Mki67"), pt.size = 0.1) 
FeaturePlot(ENS_Heart_051923.QC.1, features = c('Ednrb','Phox2b'),label = TRUE)
ENS_Heart_051923.QC.1.Neurons <- subset(x = ENS_Heart_051923.QC.1, ident = c(2,5,6,15))
DimPlot(object = ENS_Heart_051923.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Heart_051923.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Sox10'),label = TRUE)
FeaturePlot(ENS_Heart_051923.QC.1.Neurons, features = c('Hoxa5','Phox2b','Mki67','Ednrb'),label = TRUE)
VlnPlot(object = ENS_Heart_051923.QC.1.Neurons, features = c("Snap25",'Phox2b',"Ednrb"), pt.size = 0.1) 
FeaturePlot(ENS_Heart_051923.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Sox10'),label = F)
ENS_Heart_051923.QC.1.Phox2b <- WhichCells(object = ENS_Heart_051923.QC.1.Neurons)
ENS_Heart_051923.QC.Phox2b <- subset(ENS_Heart_051923.QC, cells = c(ENS_Heart_051923.QC.1.Phox2b))
DimPlot(object = ENS_Heart_051923.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Heart_051923.QC.Phox2b 
saveRDS(ENS_Heart_051923.QC.Phox2b, file = "ENS_Heart_051923.QC.Phox2b.rds")

ENS_Heart_090523.QC.1 <- NormalizeData(object = ENS_Heart_090523.QC, verbose = FALSE)
ENS_Heart_090523.QC.1 <- FindVariableFeatures(object = ENS_Heart_090523.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Heart_090523.QC.1))
ENS_Heart_090523.QC.1 <- ScaleData(object = ENS_Heart_090523.QC.1, verbose = FALSE)
ENS_Heart_090523.QC.1 <- RunPCA(object = ENS_Heart_090523.QC.1, verbose = FALSE)
ENS_Heart_090523.QC.1 <- RunUMAP(object = ENS_Heart_090523.QC.1, dims = 1:50)
ENS_Heart_090523.QC.1 <- FindNeighbors(object = ENS_Heart_090523.QC.1, dims = 1:50)
ENS_Heart_090523.QC.1 <- FindClusters(object = ENS_Heart_090523.QC.1, resolution = 0.6)
DefaultAssay(object = ENS_Heart_090523.QC.1) <- "RNA"
VlnPlot(object = ENS_Heart_090523.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 4,7,13,16 are Phox2b+
VlnPlot(object = ENS_Heart_090523.QC.1, features = c("Sox10"), pt.size = 0.1) #16
VlnPlot(object = ENS_Heart_090523.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 4,7,13 are Snap25+
VlnPlot(object = ENS_Heart_090523.QC.1, features = c("Mki67"), pt.size = 0.1) 
FeaturePlot(ENS_Heart_090523.QC.1, features = c('Mki67','Phox2b','Ednrb','Sox10'),label = F)
DimPlot(object = ENS_Heart_090523.QC.1, reduction = 'umap', label = TRUE)
ENS_Heart_090523.QC.1.Neurons <- subset(x = ENS_Heart_090523.QC.1, ident = c(4,7,13,16))
DimPlot(object = ENS_Heart_090523.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Heart_090523.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Sox10'),label = TRUE)
FeaturePlot(ENS_Heart_090523.QC.1.Neurons, features = c('Hoxa5','Gsta4','Tbx3','Gata3'),label = TRUE)
VlnPlot(object = ENS_Heart_090523.QC.1.Neurons, features = c("Snap25",'Phox2b'), pt.size = 0.1) 
ENS_Heart_090523.QC.1.Phox2b <- WhichCells(object = ENS_Heart_090523.QC.1.Neurons)
ENS_Heart_090523.QC.Phox2b <- subset(ENS_Heart_090523.QC, cells = c(ENS_Heart_090523.QC.1.Phox2b))
DimPlot(object = ENS_Heart_090523.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Heart_090523.QC.Phox2b 
saveRDS(ENS_Heart_090523.QC.Phox2b, file = "ENS_Heart_090523.QC.Phox2b.rds")

#ENS+Gut
ENS_Gut_090523 <- Read10X(data.dir = "ENS_Gut_090523")
ENS_Gut_090523 <- CreateSeuratObject(counts = ENS_Gut_090523, project = "ENS_Gut_090523")
ENS_Gut_090523 #2412

ENS_Gut_081423 <- Read10X(data.dir = "ENS_Gut_081423")
ENS_Gut_081423 <- CreateSeuratObject(counts = ENS_Gut_081423, project = "ENS_Gut_081423")
ENS_Gut_081423 #1956

ENS_Gut_051923 <- Read10X(data.dir = "ENS_Gut_051923")
ENS_Gut_051923 <- CreateSeuratObject(counts = ENS_Gut_051923, project = "ENS_Gut_051923")
ENS_Gut_051923 #2613 cells

ENS_Gut_042023 <- Read10X(data.dir = "ENS_Gut_042023")
ENS_Gut_042023 <- CreateSeuratObject(counts = ENS_Gut_042023, project = "ENS_Gut_042023")
ENS_Gut_042023 #6405

## Quanlity control 
ENS_Gut_090523
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Gut_090523), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Gut_090523, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Gut_090523, slot = 'counts'))
ENS_Gut_090523[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Gut_090523, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_090523.QC <- subset(x = ENS_Gut_090523, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Gut_090523.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_090523.QC #2369

ENS_Gut_081423
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Gut_081423), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Gut_081423, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Gut_081423, slot = 'counts'))
ENS_Gut_081423[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Gut_081423, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_081423.QC <- subset(x = ENS_Gut_081423, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Gut_081423.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_081423.QC #1928

ENS_Gut_051923
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Gut_051923), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Gut_051923, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Gut_051923, slot = 'counts'))
ENS_Gut_051923[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Gut_051923, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_051923.QC <- subset(x = ENS_Gut_051923, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Gut_051923.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_051923.QC #2564

ENS_Gut_042023
mito.features <- grep(pattern = "^mt-", x = rownames(x = ENS_Gut_042023), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ENS_Gut_042023, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ENS_Gut_042023, slot = 'counts'))
ENS_Gut_042023[['percent.mito']] <- percent.mito
VlnPlot(object = ENS_Gut_042023, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_042023.QC <- subset(x = ENS_Gut_042023, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = ENS_Gut_042023.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
ENS_Gut_042023.QC #6383

## rename datasets, add sample id
ENS_Gut_090523.QC <- RenameCells(ENS_Gut_090523.QC, add.cell.id = "ENS_Gut_090523")
Cells(ENS_Gut_090523.QC)
ENS_Gut_090523.QC

ENS_Gut_081423.QC <- RenameCells(ENS_Gut_081423.QC, add.cell.id = "ENS_Gut_081423")
Cells(ENS_Gut_081423.QC)
ENS_Gut_081423.QC

ENS_Gut_051923.QC <- RenameCells(ENS_Gut_051923.QC, add.cell.id = "ENS_Gut_051923")
Cells(ENS_Gut_051923.QC)
ENS_Gut_051923.QC

ENS_Gut_042023.QC <- RenameCells(ENS_Gut_042023.QC, add.cell.id = "ENS_Gut_042023")
Cells(ENS_Gut_042023.QC)
ENS_Gut_042023.QC

#Unsupervised clustering
ENS_Gut_090523.QC.1 <- NormalizeData(object = ENS_Gut_090523.QC, verbose = FALSE)
ENS_Gut_090523.QC.1 <- FindVariableFeatures(object = ENS_Gut_090523.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_090523.QC.1))
ENS_Gut_090523.QC.1 <- ScaleData(object = ENS_Gut_090523.QC.1, verbose = FALSE)
ENS_Gut_090523.QC.1 <- RunPCA(object = ENS_Gut_090523.QC.1, verbose = FALSE)
ENS_Gut_090523.QC.1 <- RunUMAP(object = ENS_Gut_090523.QC.1, dims = 1:50)
ENS_Gut_090523.QC.1 <- FindNeighbors(object = ENS_Gut_090523.QC.1, dims = 1:50)
ENS_Gut_090523.QC.1 <- FindClusters(object = ENS_Gut_090523.QC.1, resolution = 0.6)
DefaultAssay(object = ENS_Gut_090523.QC.1) <- "RNA"
VlnPlot(object = ENS_Gut_090523.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 5 are Phox2b+
VlnPlot(object = ENS_Gut_090523.QC.1, features = c("Sox10"), pt.size = 0.1) #None
VlnPlot(object = ENS_Gut_090523.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 5 are Snap25+
VlnPlot(object = ENS_Gut_090523.QC.1, features = c("Mki67"), pt.size = 0.1)  
FeaturePlot(ENS_Gut_090523.QC.1, features = c('Acta2','Phox2b',"Myh11","Epcam"),label = TRUE)
ENS_Gut_090523.QC.1.Neurons <- subset(x = ENS_Gut_090523.QC.1, ident = c(5))
DimPlot(object = ENS_Gut_090523.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Gut_090523.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Sox10'),label = TRUE)
FeaturePlot(ENS_Gut_090523.QC.1.Neurons, features = c('Hoxa5','Phox2b','Mki67','Snap25'),label = TRUE)
VlnPlot(object = ENS_Gut_090523.QC.1.Neurons, features = c("Snap25",'Phox2b'), pt.size = 0.1) 
ENS_Gut_090523.QC.1.Phox2b <- WhichCells(object = ENS_Gut_090523.QC.1.Neurons)
ENS_Gut_090523.QC.Phox2b <- subset(ENS_Gut_090523.QC, cells = c(ENS_Gut_090523.QC.1.Phox2b))
DimPlot(object = ENS_Gut_090523.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Gut_090523.QC.Phox2b #205 cells
saveRDS(ENS_Gut_090523.QC.Phox2b, file = "ENS_Gut_090523.QC.Phox2b.rds")

ENS_Gut_081423.QC.1 <- NormalizeData(object = ENS_Gut_081423.QC, verbose = FALSE)
ENS_Gut_081423.QC.1 <- FindVariableFeatures(object = ENS_Gut_081423.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_081423.QC.1))
ENS_Gut_081423.QC.1 <- ScaleData(object = ENS_Gut_081423.QC.1, verbose = FALSE)
ENS_Gut_081423.QC.1 <- RunPCA(object = ENS_Gut_081423.QC.1, verbose = FALSE)
ENS_Gut_081423.QC.1 <- RunUMAP(object = ENS_Gut_081423.QC.1, dims = 1:50)
ENS_Gut_081423.QC.1 <- FindNeighbors(object = ENS_Gut_081423.QC.1, dims = 1:50)
ENS_Gut_081423.QC.1 <- FindClusters(object = ENS_Gut_081423.QC.1, resolution = 0.6)
DefaultAssay(object = ENS_Gut_081423.QC.1) <- "RNA"
VlnPlot(object = ENS_Gut_081423.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 1,4,5,6,9 are Phox2b+
VlnPlot(object = ENS_Gut_081423.QC.1, features = c("Sox10"), pt.size = 0.1) #None
VlnPlot(object = ENS_Gut_081423.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 1,4,5,6 are Snap25+
VlnPlot(object = ENS_Gut_081423.QC.1, features = c("Mki67"), pt.size = 0.1)  
FeaturePlot(ENS_Gut_081423.QC.1, features = c('Ednrb','Phox2b',"Snap25","Mki67"),label = TRUE)
ENS_Gut_081423.QC.1.Neurons <- subset(x = ENS_Gut_081423.QC.1, ident = c(1,4,5,6,9)) #c9 needs to be checked
DimPlot(object = ENS_Gut_081423.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Gut_081423.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Sox10'),label = TRUE) # some cells are phox2b negative
FeaturePlot(ENS_Gut_081423.QC.1.Neurons, features = c('Ednrb','Phox2b','Mki67','Snap25'),label = TRUE)
VlnPlot(object = ENS_Gut_081423.QC.1.Neurons, features = c("Snap25",'Phox2b'), pt.size = 0.1) 
ENS_Gut_081423.QC.1.Phox2b <- WhichCells(object = ENS_Gut_081423.QC.1.Neurons) 
ENS_Gut_081423.QC.Phox2b <- subset(ENS_Gut_081423.QC, cells = c(ENS_Gut_081423.QC.1.Phox2b))
DimPlot(object = ENS_Gut_081423.QC.Phox2b, reduction = 'umap', label = TRUE)

ENS_Gut_081423.QC.Phox2b.1 <- NormalizeData(object = ENS_Gut_081423.QC.Phox2b, verbose = FALSE)
ENS_Gut_081423.QC.Phox2b.1 <- FindVariableFeatures(object = ENS_Gut_081423.QC.Phox2b.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_081423.QC.Phox2b.1))
ENS_Gut_081423.QC.Phox2b.1 <- ScaleData(object = ENS_Gut_081423.QC.Phox2b.1, verbose = FALSE)
ENS_Gut_081423.QC.Phox2b.1 <- RunPCA(object = ENS_Gut_081423.QC.Phox2b.1, verbose = FALSE)
ENS_Gut_081423.QC.Phox2b.1 <- RunUMAP(object = ENS_Gut_081423.QC.Phox2b.1, dims = 1:50)
ENS_Gut_081423.QC.Phox2b.1 <- FindNeighbors(object = ENS_Gut_081423.QC.Phox2b.1, dims = 1:50)
ENS_Gut_081423.QC.Phox2b.1 <- FindClusters(object = ENS_Gut_081423.QC.Phox2b.1, resolution = 4)
ENS_Gut_081423.QC.Phox2b.1 
DimPlot(ENS_Gut_081423.QC.Phox2b.1, label = T)
VlnPlot(object = ENS_Gut_081423.QC.Phox2b.1, features = c("Snap25",'Phox2b'), pt.size = 0.1) #c3? c16? -> Phox2b low -> remove
FeaturePlot(ENS_Gut_081423.QC.Phox2b.1, features = c('Ednrb','Phox2b','Mki67','Snap25'),label = TRUE)
ENS_Gut_081423.QC.Phox2b.1 <- subset(x = ENS_Gut_081423.QC.Phox2b.1, ident = c(0,1,2,4,5,6,7,8,9,10,11,12,13,14,15,17,18)) 
DimPlot(ENS_Gut_081423.QC.Phox2b.1, label = T)
ENS_Gut_081423.QC.1.Phox2b <- WhichCells(object = ENS_Gut_081423.QC.Phox2b.1) 
ENS_Gut_081423.QC.Phox2b <- subset(ENS_Gut_081423.QC, cells = c(ENS_Gut_081423.QC.1.Phox2b)) 
DimPlot(ENS_Gut_081423.QC.Phox2b)
saveRDS(ENS_Gut_081423.QC.Phox2b, "ENS_Gut_081423.QC.Phox2b.rds")

ENS_Gut_051923.QC.1 <- NormalizeData(object = ENS_Gut_051923.QC, verbose = FALSE)
ENS_Gut_051923.QC.1 <- FindVariableFeatures(object = ENS_Gut_051923.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_051923.QC.1))
ENS_Gut_051923.QC.1 <- ScaleData(object = ENS_Gut_051923.QC.1, verbose = FALSE)
ENS_Gut_051923.QC.1 <- RunPCA(object = ENS_Gut_051923.QC.1, verbose = FALSE)
ENS_Gut_051923.QC.1 <- RunUMAP(object = ENS_Gut_051923.QC.1, dims = 1:50)
ENS_Gut_051923.QC.1 <- FindNeighbors(object = ENS_Gut_051923.QC.1, dims = 1:50)
ENS_Gut_051923.QC.1 <- FindClusters(object = ENS_Gut_051923.QC.1, resolution = 4) #0.6 -> 4
DefaultAssay(object = ENS_Gut_051923.QC.1) <- "RNA"
VlnPlot(object = ENS_Gut_051923.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 1,3,7,8,9,10? are Phox2b+
#resolution = 4: 0,3,6,7,8,17,19,23
VlnPlot(object = ENS_Gut_051923.QC.1, features = c("Sox10"), pt.size = 0.1) #None
VlnPlot(object = ENS_Gut_051923.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 1,3,7,8,9 are Snap25+
VlnPlot(object = ENS_Gut_051923.QC.1, features = c("Mki67"), pt.size = 0.1) #5,7 
FeaturePlot(ENS_Gut_051923.QC.1, features = c('Acta2','Phox2b',"Myh11","Mylk","Ednrb","Mki67"),label = TRUE) 
ENS_Gut_051923.QC.1.Neurons <- subset(x = ENS_Gut_051923.QC.1, ident = c(0,3,6,7,8,17,19,23))
DimPlot(object = ENS_Gut_051923.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Gut_051923.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Ednrb'),label = TRUE)
FeaturePlot(ENS_Gut_051923.QC.1.Neurons, features = c('Sox10','Phox2b','Mki67','Mylk'),label = TRUE)
VlnPlot(object = ENS_Gut_051923.QC.1.Neurons, features = c("Snap25",'Phox2b'), pt.size = 0.1) 
ENS_Gut_051923.QC.1.Phox2b <- WhichCells(object = ENS_Gut_051923.QC.1.Neurons)
ENS_Gut_051923.QC.Phox2b <- subset(ENS_Gut_051923.QC, cells = c(ENS_Gut_051923.QC.1.Phox2b))
DimPlot(object = ENS_Gut_051923.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Gut_051923.QC.Phox2b 
# cluster again with high resolution to remove phox2b negative cells
ENS_Gut_051923.QC.Phox2b.1 <- NormalizeData(object = ENS_Gut_051923.QC.Phox2b, verbose = FALSE)
ENS_Gut_051923.QC.Phox2b.1 <- FindVariableFeatures(object = ENS_Gut_051923.QC.Phox2b.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_051923.QC.Phox2b.1))
ENS_Gut_051923.QC.Phox2b.1 <- ScaleData(object = ENS_Gut_051923.QC.Phox2b.1, verbose = FALSE)
ENS_Gut_051923.QC.Phox2b.1 <- RunPCA(object = ENS_Gut_051923.QC.Phox2b.1, verbose = FALSE)
ENS_Gut_051923.QC.Phox2b.1 <- RunUMAP(object = ENS_Gut_051923.QC.Phox2b.1, dims = 1:50)
ENS_Gut_051923.QC.Phox2b.1 <- FindNeighbors(object = ENS_Gut_051923.QC.Phox2b.1, dims = 1:50)
ENS_Gut_051923.QC.Phox2b.1 <- FindClusters(object = ENS_Gut_051923.QC.Phox2b.1, resolution = 4)
DimPlot(ENS_Gut_051923.QC.Phox2b.1, label = T)
FeaturePlot(ENS_Gut_051923.QC.Phox2b.1, features = c('Sox10','Phox2b','Mki67','Ednrb'),label = TRUE)
VlnPlot(ENS_Gut_051923.QC.Phox2b.1, features = c('Sox10','Phox2b','Mki67','Ednrb'))  
ENS_Gut_051923.QC.1.Neurons <- ENS_Gut_051923.QC.Phox2b.1
DimPlot(ENS_Gut_051923.QC.1.Neurons)
FeaturePlot(ENS_Gut_051923.QC.1.Neurons, features = c('Sox10','Phox2b','Mki67','Ednrb'),label = TRUE)
ENS_Gut_051923.QC.1.Phox2b <- WhichCells(object = ENS_Gut_051923.QC.1.Neurons)
ENS_Gut_051923.QC.Phox2b <- subset(ENS_Gut_051923.QC, cells = c(ENS_Gut_051923.QC.1.Phox2b))
DimPlot(object = ENS_Gut_051923.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Gut_051923.QC.Phox2b 
saveRDS(ENS_Gut_051923.QC.Phox2b, file = "ENS_Gut_051923.QC.Phox2b.rds")

ENS_Gut_042023.QC.1 <- NormalizeData(object = ENS_Gut_042023.QC, verbose = FALSE)
ENS_Gut_042023.QC.1 <- FindVariableFeatures(object = ENS_Gut_042023.QC.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_042023.QC.1))
ENS_Gut_042023.QC.1 <- ScaleData(object = ENS_Gut_042023.QC.1, verbose = FALSE)
ENS_Gut_042023.QC.1 <- RunPCA(object = ENS_Gut_042023.QC.1, verbose = FALSE)
ENS_Gut_042023.QC.1 <- RunUMAP(object = ENS_Gut_042023.QC.1, dims = 1:50)
ENS_Gut_042023.QC.1 <- FindNeighbors(object = ENS_Gut_042023.QC.1, dims = 1:50)
ENS_Gut_042023.QC.1 <- FindClusters(object = ENS_Gut_042023.QC.1, resolution = 10) #0.6 to 4 to 6 to 10
DefaultAssay(object = ENS_Gut_042023.QC.1) <- "RNA"
VlnPlot(object = ENS_Gut_042023.QC.1, features = c("Phox2b"), pt.size = 0.1) #Cluster 1,2,4,10 are Phox2b+
#resolution = 4: 1,2,10?,11,13,14,16,17,18,19,22,24,26,28,30
#resolution = 6: 0,3,4,5?,11,12,16,17,18,22,24,25,29,30,32,33,35,36
#resolution = 10: 0,1?,2,4,6,8,9,10,11,16,17,18,19,20,32,36,38,39,41,44,45,47,52?,53?,55,58
VlnPlot(object = ENS_Gut_042023.QC.1, features = c("Sox10"), pt.size = 0.1) #None
VlnPlot(object = ENS_Gut_042023.QC.1, features = c("Snap25"), pt.size = 0.1) #Cluster 1,2,4,10 are Snap25+
VlnPlot(object = ENS_Gut_042023.QC.1, features = c("Mki67"), pt.size = 0.1)  
FeaturePlot(ENS_Gut_042023.QC.1, features = c('Phox2b','Ednrb'),label = TRUE)
ENS_Gut_042023.QC.1.Neurons <- subset(x = ENS_Gut_042023.QC.1, ident = c( 0,2,4,6,8,9,10,11,16,17,18,19,20,32,36,38,39,41,44,45,47,55,58))
DimPlot(object = ENS_Gut_042023.QC.1.Neurons, reduction = 'umap', label = TRUE)
FeaturePlot(ENS_Gut_042023.QC.1.Neurons, features = c('Phox2b','Ednrb'),label = TRUE)
FeaturePlot(ENS_Gut_042023.QC.1.Neurons, features = c('Snap25','Phox2b','Mki67','Ednrb'),label = TRUE)
VlnPlot(object = ENS_Gut_042023.QC.1.Neurons, features = c("Phox2b"), pt.size = 0.1)  
FeaturePlot(ENS_Gut_042023.QC.1.Neurons, features = c("Actg2",'Phox2b','Myh11','Mylk'),label = TRUE) 
ENS_Gut_042023.QC.1.Phox2b <- WhichCells(object = ENS_Gut_042023.QC.1.Neurons) 
ENS_Gut_042023.QC.Phox2b <- subset(ENS_Gut_042023.QC, cells = c(ENS_Gut_042023.QC.1.Phox2b))
DimPlot(object = ENS_Gut_042023.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Gut_042023.QC.Phox2b 
ENS_Gut_042023.QC.Phox2b.1 <- NormalizeData(object = ENS_Gut_042023.QC.Phox2b, verbose = FALSE)
ENS_Gut_042023.QC.Phox2b.1 <- FindVariableFeatures(object = ENS_Gut_042023.QC.Phox2b.1, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = ENS_Gut_042023.QC.Phox2b.1))
ENS_Gut_042023.QC.Phox2b.1 <- ScaleData(object = ENS_Gut_042023.QC.Phox2b.1, verbose = FALSE)
ENS_Gut_042023.QC.Phox2b.1 <- RunPCA(object = ENS_Gut_042023.QC.Phox2b.1, verbose = FALSE) 
ENS_Gut_042023.QC.Phox2b.1 <- RunUMAP(object = ENS_Gut_042023.QC.Phox2b.1, dims = 1:50)
ENS_Gut_042023.QC.Phox2b.1 <- FindNeighbors(object = ENS_Gut_042023.QC.Phox2b.1, dims = 1:50)
ENS_Gut_042023.QC.Phox2b.1 <- FindClusters(object = ENS_Gut_042023.QC.Phox2b.1, resolution = 4)
DimPlot(ENS_Gut_042023.QC.Phox2b.1)
FeaturePlot(ENS_Gut_042023.QC.Phox2b.1, features = c('Actg2'),label = TRUE)
FeaturePlot(ENS_Gut_042023.QC.Phox2b.1, features = c("nCount_RNA","nFeature_RNA","percent.mito"),label = TRUE)
VlnPlot(ENS_Gut_042023.QC.Phox2b.1, features = c('Ednrb','Phox2b',"Mki67","Sox10"))
FeaturePlot(ENS_Gut_042023.QC.Phox2b.1, features = c('Ednrb','Phox2b',"Col1a2","Itga1","Snap25"))
ENS_Gut_042023.QC.1.Phox2b <- WhichCells(object = ENS_Gut_042023.QC.Phox2b.1)
ENS_Gut_042023.QC.Phox2b <- subset(ENS_Gut_042023.QC, cells = c(ENS_Gut_042023.QC.1.Phox2b))
DimPlot(object = ENS_Gut_042023.QC.Phox2b, reduction = 'umap', label = TRUE)
ENS_Gut_042023.QC.Phox2b 
saveRDS(ENS_Gut_042023.QC.Phox2b, "ENS_Gut_042023.QC.Phox2b.rds")

ENS_Heart_051923 <- readRDS("ENS_Heart_051923.QC.Phox2b.rds") 
ENS_Heart_090523 <- readRDS("ENS_Heart_090523.QC.Phox2b.rds") 
ENS_Gut_051923 <- readRDS("ENS_Gut_051923.QC.Phox2b.rds") 
ENS_Gut_081423 <- readRDS("ENS_Gut_081423.QC.Phox2b.rds") 
ENS_Gut_090523 <- readRDS("ENS_Gut_090523.QC.Phox2b.rds") 
ENS_Gut_042023 <- readRDS("ENS_Gut_042023.QC.Phox2b.rds") 

## Apply more stringent QC for multiple dataset integration





            



            
