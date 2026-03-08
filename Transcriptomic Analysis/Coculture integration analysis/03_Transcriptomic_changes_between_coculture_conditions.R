</> R
## R version 4.1.2 (2021-11-01)
## SeuratObject_4.1.3 sp_2.2-0    

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

## Compare ICNS vs ENS neurons to define system-specific marker genes. Use that to capture changes induced by ENS-heart

#ENS_ICNS.rds (ENS_ICNS) is the integration of the E16.5, E18.5 ICNS scRNA-seq dataset of this study with 
#the E15.5 and E18.5 ENS scRNA-seq dataset from GSE149524, Morarach, K., et al., Nat Neurosci, 2021.
# "ICN.neuron" and "ENS.neuron" include only the neuron populations.
# "ICN.precursors" and "ENS.precursors" include only the precursor populations.
# "ICN.neuroblast" and "ENS.neuroblast" include only the neuroblast populations. 
# "common_genes" includes genes that are commonly expressed by our scRNA-seq datasets and GSE149524 datasets, since a different
# version of CellRanger processed GSE149524 datasets.

# Acquire DEGs in the precursor or neuroblast state using E15.5 ENS and E16.5 ICNS
sub <- subset(ENS_ICNS, orig.ident %in% c("EN.E15","ICN.E16.2"))
DefaultAssay(sub) <- "RNA"
neuroblast.d <- FindMarkers(sub, ident.1 = 'ICN.neuroblast', ident.2 = 'ENS.neuroblast', only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
neuroblast.d$DE_group <- ifelse(neuroblast.d$avg_log2FC > 0,
                            "ICN",
                            "ENS")
neuroblast.d <- tibble::rownames_to_column(neuroblast.d, var = "gene")
write.table(neuroblast.d, "ENS_E15_ICN_E16_neuroblasts.xls",  sep = "\t")

precursor.d <- FindMarkers(sub, ident.1 = 'ICN.precursors', ident.2 = "ENS.precursors", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
precursor.d$DE_group <- ifelse(precursor.d$avg_log2FC > 0,
                             "ICN",
                             "ENS")
precursor.d <- tibble::rownames_to_column(precursor.d, var = "gene")
write.table(precursor.d, "ENS_E15_ICN_E16_precurosrs.xls",  sep = "\t")

# Acquire DEGs in the precursor or neuroblast state using E18.5 ENS and E18.5 ICNS
neuron.d <- FindMarkers(ENS_ICNS_E18.5, ident.1 = 'ICN.neuron', ident.2 = "ENS.neuron", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
neuron.d$DE_group <- ifelse(neuron.d$avg_log2FC > 0,
                             "ICN",
                             "ENS")
neuron.d <- tibble::rownames_to_column(neuron.d, var = "gene")
write.table(neuron.d, "ENS_E18_ICN_E18_neuron.xls",  sep = "\t")


## Compare neurons of ENS-heart and ENS-intestine for DEGs 
n <- subset(coculture, state %in% "Neurons")
DefaultAssay(n) <- "RNA"
d <- FindMarkers(n, ident.1 = 'ENS_Heart', ident.2 = 'ENS_Gut', only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
d$DE_group <- ifelse(d$avg_log2FC > 0,
                     "ENS_Heart",
                     "ENS_Gut")
d <- tibble::rownames_to_column(d, var = "gene")
write.table(d, "coculture_Neurons_deg.xls",  sep = "\t")

neuron_heart <- neuron.d %>% filter(DE_group == "ICN",
                                  p_val_adj < 0.05) 
neuron_gut <- neuron.d %>% filter(DE_group == "ENS",
                                p_val_adj < 0.05) 
cc_neuron_heart <- d %>% filter(DE_group == "ENS_Heart",
                                        p_val_adj < 0.05) 
cc_neuron_gut <- d %>% filter(DE_group == "ENS_Gut",
                                      p_val_adj < 0.05) 



############### Compute ICNS-ENS similarity scores ############### 

neuron_heart <- neuron %>% filter(DE_group == "ICN",
                                  p_val_adj < 0.05) #1636
neuron_gut <- neuron %>% filter(DE_group == "ENS",
                                p_val_adj < 0.05) #833

neuroblast_heart <- neuroblast %>% filter(DE_group == "ICN",
                                          p_val_adj < 0.05) #1492
neuroblast_gut <- neuroblast %>% filter(DE_group == "ENS",
                                        p_val_adj < 0.05) #816

precursor_heart <- precursor %>% filter(DE_group == "ICN",
                                        p_val_adj < 0.05) #2135
precursor_gut <- precursor %>% filter(DE_group == "ENS",
                                      p_val_adj < 0.05) #943


#### Prep the gene and cells to test ####
Idents(cc) <- "state"
DimPlot(cc)
sub <- subset(cc, state %in% c("Neuroblast"))
DimPlot(sub, split.by = "ID")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- neuroblast_gut$gene

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_neuroblast",
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:825]
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_neuroblast"]] <- meta[,1]
FeaturePlot(sub,"ens_neuroblast", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_neuroblast - sub$ens_neuroblast

FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "neuroblast_score.rds")

Idents(cc) <- "state"
DimPlot(cc)
sub <- subset(cc, state %in% c("Neurons"))
DimPlot(sub, split.by = "ID")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- neuron_gut$gene

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_neuron",
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:842]
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_neuron"]] <- meta[,1]
FeaturePlot(sub,"ens_neuron", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_neuron - sub$ens_neuron

FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "neuron_score.rds")

Idents(cc) <- "state"
DimPlot(cc)
sub <- subset(cc, state %in% c("Precursors"))
DimPlot(sub, split.by = "ID")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- precursor_gut$gene

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_precursor",
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:952]
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_precursor"]] <- meta[,1]
FeaturePlot(sub,"ens_precursor", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_precursor - sub$ens_precursor
FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "precursor_score.rds")




