rm(list = ls())
library(Seurat)
library(SeuratDisk)

# Neuron cells
data_path <- "/gpfs/gibbs/project/zhao/jz874/rui-lab/Yvonne-ICNs/from_Yvonne/2024.03.11/Organ_neurons"
oin <-  readRDS(paste0(data_path, "/E14.5_OIN.rds"))
DefaultAssay(oin) <- "RNA"

# heart
heart_neuron <- oin[, as.vector(oin@meta.data$orig.ident) == "ICN.E14.3"]
# gut
gut_neuron <- oin[, as.vector(oin@meta.data$orig.ident) == "ENS_E14"]
# pancrease
pan_neuron <- oin[, as.vector(oin@meta.data$orig.ident) == "E14.5_pancrease"]
# lung
lung_neuron <- oin[, as.vector(oin@meta.data$orig.ident) == "E14.5_Lung"]

# create new objs
heart_obj <- CreateSeuratObject(counts = heart_neuron@assays$RNA@counts, meta.data = heart_neuron@meta.data)
rm(heart_neuron)
gut_obj <- CreateSeuratObject(counts = gut_neuron@assays$RNA@counts, meta.data = gut_neuron@meta.data)
rm(gut_neuron)
pan_obj <- CreateSeuratObject(counts = pan_neuron@assays$RNA@counts, meta.data = pan_neuron@meta.data)
rm(pan_neuron)
lung_obj <- CreateSeuratObject(counts = lung_neuron@assays$RNA@counts, meta.data = lung_neuron@meta.data)
rm(lung_neuron)

# save data
setwd("/gpfs/gibbs/project/zhao/jz874/rui-lab/Yvonne-ICNs/oin_trajectory_integration_portal/preprocessed_data")
SaveH5Seurat(heart_obj, filename = "oin_heart.h5Seurat")
Convert("oin_heart.h5Seurat", dest = "h5ad")
SaveH5Seurat(gut_obj, filename = "oin_gut.h5Seurat")
Convert("oin_gut.h5Seurat", dest = "h5ad")
SaveH5Seurat(pan_obj, filename = "oin_pan.h5Seurat")
Convert("oin_pan.h5Seurat", dest = "h5ad")
SaveH5Seurat(lung_obj, filename = "oin_lung.h5Seurat")
Convert("oin_lung.h5Seurat", dest = "h5ad")

