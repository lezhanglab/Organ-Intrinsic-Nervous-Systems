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
# "common_genes" includes genes that are commonly present in our scRNA-seq datasets and GSE149524 datasets, since a different
# version of CellRanger processed GSE149524 datasets.

# 1. Acquire DEGs in the precursor or neuroblast state using E15.5 ENS and E16.5 ICNS
sub <- subset(ENS_ICNS, orig.ident %in% c("ENS.E15.5","ICNS.E16.5"))
DefaultAssay(sub) <- "RNA"
neuroblast_d <- FindMarkers(sub, ident.1 = 'ICN.neuroblast', ident.2 = 'ENS.neuroblast', only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
neuroblast_d$DE_group <- ifelse(neuroblast_d$avg_log2FC > 0,
                            "ICN",
                            "ENS")
neuroblast_d <- tibble::rownames_to_column(neuroblast_d, var = "gene")
write.table(neuroblast_d, "ENS_E15_ICN_E16_neuroblasts.xls",  sep = "\t")

precursor_d <- FindMarkers(sub, ident.1 = 'ICN.precursors', ident.2 = "ENS.precursors", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
precursor_d$DE_group <- ifelse(precursor_d$avg_log2FC > 0,
                             "ICN",
                             "ENS")
precursor_d <- tibble::rownames_to_column(precursor_d, var = "gene")
write.table(precursor_d, "ENS_E15_ICN_E16_precursors.xls",  sep = "\t")

# Acquire DEGs in the precursor or neuroblast state using E18.5 ENS and E18.5 ICNS
sub <- subset(ENS_ICNS, orig.ident %in% c("ENS.E18.5","ICNS.E18.5"))
neuron_d <- FindMarkers(sub, ident.1 = 'ICN.neuron', ident.2 = "ENS.neuron", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
neuron_d$DE_group <- ifelse(neuron_d$avg_log2FC > 0,
                             "ICN",
                             "ENS")
neuron_d <- tibble::rownames_to_column(neuron_d, var = "gene")
write.table(neuron_d, "ENS_E18_ICN_E18_neuron.xls",  sep = "\t")


## Compare neurons of ENS-heart and ENS-intestine for DEGs 
n <- subset(coculture, state %in% "Neurons")
DefaultAssay(n) <- "RNA"
d <- FindMarkers(n, ident.1 = 'ENS_Heart', ident.2 = 'ENS_Gut', only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25, test.use = 'wilcox', features = common_genes)
d$DE_group <- ifelse(d$avg_log2FC > 0,
                     "ENS_Heart",
                     "ENS_Gut")
d <- tibble::rownames_to_column(d, var = "gene")
write.table(d, "coculture_Neurons_deg.xls",  sep = "\t")

neuron_heart <- neuron_d %>% filter(DE_group == "ICN",
                                  p_val_adj < 0.05) 
neuron_gut <- neuron_d %>% filter(DE_group == "ENS",
                                p_val_adj < 0.05) 
cc_neuron_heart <- d %>% filter(DE_group == "ENS_Heart",
                                        p_val_adj < 0.05) 
cc_neuron_gut <- d %>% filter(DE_group == "ENS_Gut",
                                      p_val_adj < 0.05) 
# Use neuron_heart and cc_neuron_heart, we found that of the 663 genes upregulated in ENS-heart relative to ENS-intestine at the neuron state, 
# 372 (56%) were also upregulated in ICNS relative to ENS neurons at E18.5. 

#2. Compute ICNS-ENS similarity scores at the precursor, neuronblast or neuron state.

neuron_heart <- neuron_d %>% filter(DE_group == "ICN",
                                  p_val_adj < 0.05) 
neuron_gut <- neuron_d %>% filter(DE_group == "ENS",
                                p_val_adj < 0.05) 

neuroblast_heart <- neuroblast_d %>% filter(DE_group == "ICN",
                                          p_val_adj < 0.05) 
neuroblast_gut <- neuroblast_d %>% filter(DE_group == "ENS",
                                        p_val_adj < 0.05) 

precursor_heart <- precursor_d %>% filter(DE_group == "ICN",
                                        p_val_adj < 0.05) 
precursor_gut <- precursor_d %>% filter(DE_group == "ENS",
                                      p_val_adj < 0.05) 


#### Prepare the gene and cells to test ####
Idents(coculture) <- "state"
sub <- subset(coculture, state %in% c("Neuroblast"))
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- neuroblast_gut$gene #change this to "neuroblast_heart" for "icn_neuroblast"

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_neuroblast", #change this to "icn_neuroblast" for "icn_neuroblast"
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:825] #change accroding to gene number 
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_neuroblast"]] <- meta[,1]
FeaturePlot(sub,"ens_neuroblast", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_neuroblast - sub$ens_neuroblast

FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "neuroblast_score.rds")


sub <- subset(coculture, state %in% c("Neurons")) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- neuron_gut$gene #change this to "neuron_heart" for "icns_neuron"

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_neuron", #change this to "icns_neuron" for "icns_neuron"
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:842] #change accroding to gene number 
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_neuron"]] <- meta[,1]
FeaturePlot(sub,"ens_neuron", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_neuron - sub$ens_neuron

FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "neuron_score.rds")


sub <- subset(coculture, state %in% c("Precursors"))
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "ID"

gene <- precursor_gut$gene #change this to "precursor_heart" for "icns_precursor"

set.seed(1)
test <- AddModuleScore(
  sub,
  features = gene,
  name = "ens_precursor", #change this to "icns_precursor" for "icns_precursor"
  ctrl = 200,
  assay = "RNA")

meta <- test@meta.data
meta <- meta[,10:952] #change accroding to gene number 
meta <- rowMeans(meta) %>% 
  as.data.frame()
sub[["ens_precursor"]] <- meta[,1]
FeaturePlot(sub,"ens_precursor", split.by = "ID")

sub$icn_ens_similarity <- sub$icn_precursor - sub$ens_precursor
FeaturePlot(sub, "icn_ens_similarity", split.by = "ID", 
            cols = rev(brewer.pal(n = 11, name = "RdBu")))
saveRDS(sub, "precursor_score.rds")




