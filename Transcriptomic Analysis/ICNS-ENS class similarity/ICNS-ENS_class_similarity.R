library(Seurat)
library(lsa)
setwd("/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/from_Yvonne/2026.01.07")
obj <- readRDS("ENS_ICNS_balance_7_annotate.rds")

DimPlot(obj, label = T)



## PCA after integration, cosine similarity
mat <- obj@reductions$pca@cell.embeddings

meta <- obj@meta.data
anno <- as.vector(meta$NeuronCluster)
anno_chosen <- c("ICNS_Npy","InM_myENC8-9","ICNS_Ddah1","IPAN_myENC6","IPAN_myENC7",
                 "myENC4","ExM_myENC1–3",
                 "IPAN_myENC12","IN_myENC10")
n_anno_type <- length(anno_chosen)  
n_dim <- dim(mat)[2]
average_anno_mat <- matrix(0, nrow = n_anno_type, ncol = n_dim)
for (i in 1:n_anno_type) {
  anno_type <- anno_chosen[i]
  print(anno_type)
  mat_anno_tmp <- mat[anno == anno_type, ]
  average_anno_mat[i, ] <- colMeans(mat_anno_tmp)
}
average_anno_mat <- t(average_anno_mat)
colnames(average_anno_mat) <- anno_chosen

similarity_mat <- cosine(average_anno_mat)

library(reshape2)
library(ggplot2)
df <- melt(similarity_mat)
p <- ggplot(df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  coord_fixed(expand = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = median(df$value, na.rm = TRUE)
  )
p
similarity_mat

setwd("/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_revision_similarity_calculation")
write.csv(similarity_mat, 'similarity_matrix_allpopulations.csv')

pdf("similarity_matrix_allpopulations.pdf", width = 7.72, height = 5.51)
p
dev.off()



