</> R
## R version 4.1.2 (2021-11-01)

library(dplyr)
library(colorRamp2)
library(ComplexHeatmap)

ss_cluster <- readRDS("ss_cluster_PrecursortoNeuroblast.rds")
gene_ts <- readRDS("PrecursortoNeuroblast_final_time_expN.rds")

# pseudotime
pt <- gene_ts$time

# expression matrix with genes as rows
data <- gene_ts[, -1] #remove the "time" column
data <- t(data) %>% as.data.frame()

# keep gene names explicitly
data$gene <- rownames(data)

# join cluster labels by gene
ss_cluster$gene <- as.character(ss_cluster$gene)
data$gene <- as.character(data$gene)
data <- left_join(data, ss_cluster, by = "gene")

# order by cluster
# define cluster order
cluster_order <- c(12,10,6,8,2,7,4,9,13,11,1,3,5)
data$ss_cluster <- factor(data$ss_cluster, levels = cluster_order)

# cluster 12,10 are grouped to cluster A
# cluster 6,8,2 are grouped to cluster B
# cluster 7,4 are grouped to cluster C
# cluster 9 are grouped to cluster D
# cluster 13 are grouped to cluster E
# cluster 11,1,3,5 are grouped to cluster F

data$label <- dplyr::case_when(
  data$ss_cluster %in% c(12,10) ~ "A",
  data$ss_cluster %in% c(6,8,2) ~ "B",
  data$ss_cluster %in% c(7,4) ~ "C",
  data$ss_cluster %in% c(9) ~ "D",
  data$ss_cluster %in% c(13) ~ "E",
  data$ss_cluster %in% c(11,1,3,5) ~ "F"
)

data$label <- factor(data$label, levels = c("A","B","C","D","E","F"))

data <- data %>% arrange(ss_cluster)

# expression matrix for plotting
mat_df <- data %>%
  select(-gene, -ss_cluster, -label)

mat <- as.matrix(mat_df)
rownames(mat) <- data$gene

genes <- c("Serpine2","Itga4","Ednrb","Sox10","Erbb3","Foxd3","Ets1","Plp1",
             "Col1a2","Dusp6","Tnc","Itgb1","Lama2","Col18a1","Ascl1",
             "Itga1","Gata3","Lef1","Tbx20","Plxna4",
             "Cux2","Sv2c","Ntrk3","Igf1r",
             "Acvr2b","Tox3","Aff2","Insm1",
             "Epha5","Cntn5","Isl1","Ncam1","Syt11",
             "L1cam","Nrp1","Notch1","Cd47","Nlgn1")
idx <- match(genes, rownames(mat))
idx
genes[is.na(idx)] #check if any missing gene

ha <- rowAnnotation(
  mark = anno_mark(
    at = idx,
    labels = genes
  )
)

# top annotation: pseudotime
pt_col_fun <- colorRamp2(
  c(0, 8.45, 16.9),
  c("#5e4fa2", "#2fb7a1", "#f1d76a")
)

colAnn <- HeatmapAnnotation(
  pseudotime = as.numeric(pt),
  col = list(pseudotime = pt_col_fun)
)

# left annotation with A-F groups, but keeps the cluster order of 12,10,6,8,... 
cluster_col <- c(
  "12" = "#E41A1C",
  "10" = "#377EB8",
  "6"  = "#4DAF4A",
  "8"  = "#984EA3",
  "2"  = "#FF7F00",
  "7"  = "#FFFF33",
  "4"  = "#A65628",
  "9"  = "#F781BF",
  "13" = "#999999",
  "11" = "#66C2A5",
  "1"  = "#FC8D62",
  "3"  = "#8DA0CB",
  "5"  = "#E78AC3"
)
rolAnn <- rowAnnotation(
  cluster = data$ss_cluster,
  col = list(cluster = cluster_col),
  show_annotation_name = FALSE
)

# heatmap colors
myCol1 <- colorRampPalette(c("steelblue1", "white", "red"))(300)
myBreaks1 <- seq(min(mat), max(mat), length.out = 300)

ht <- Heatmap(
  mat,
  col = colorRamp2(myBreaks1, myCol1),
  top_annotation = colAnn,
  left_annotation = rolAnn,
  right_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_split = data$label,
  row_title_rot = 0
)

pdf("PrecursortoNeuroblast_heatmap.pdf", width = 8, height = 10)
draw(ht)
dev.off()

