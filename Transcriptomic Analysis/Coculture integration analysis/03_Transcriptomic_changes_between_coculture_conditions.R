</> R
## R version 4.1.2 (2021-11-01)
## SeuratObject_4.1.3 sp_2.2-0    

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

## Compare ICNS vs ENS neurons at E18.5 to define system-specific marker genes
## ENS-heart vs ENS-intestine neurons to capture changes induced by heart cells. 

#ENS_ICNS_E18.5.rds is the integration of the E18.5 ICNS scRNA-seq dataset of this study with 
#the E18.5 ENS scRNA-seq dataset (GSE149524, Morarach, K., et al., Nat Neurosci, 2021).
# "ICN.neuron" and "ENS.neuron" include only the neuron populations.
neuron.d <- FindMarkers(ENS_ICNS_E18.5, ident.1 = 'ICN.neuron', ident.2 = "ENS.neuron", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'bimod')
neuron.d$DE_group <- ifelse(neuron.d$avg_log2FC > 0,
                             "ICN",
                             "ENS")
neuron.d <- tibble::rownames_to_column(neuron.d, var = "gene")
write.table(neuron.d, "ENS_E18_ICN_E18_neuron.xls",  sep = "\t")

##
DefaultAssay(n) <- "RNA"
Idents(n) <- "ID"
DimPlot(n, split.by = "ID");
d_n <- FindMarkers(n, ident.1 = 'ENS_Heart', ident.2 = "ENS_Gut", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'bimod')
d_n$DE_group <- ifelse(d_n$avg_log2FC > 0,
                        "ENS_Heart",
                        "ENS_Gut")
d_n <- tibble::rownames_to_column(d_n, var = "gene")
write.table(d_n, "coculture_d_neuron.xls",  sep = "\t")

d_n <- read.table("coculture_d_neuron.xls")



