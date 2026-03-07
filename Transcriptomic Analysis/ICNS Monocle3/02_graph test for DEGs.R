</> R
## R version 4.1.2 (2021-11-01)

library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(htmlwidgets)
library(Signac)
library(dplyr)

cds <- readRDS("PrecursortoNeuroblast.rds")
plot_cells(cds)
list_mono <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
#Check q_value cutoff
gene_list <- list_mono %>% 
  filter(q_value<0.01) # 32285 genes -> 7293 genes
saveRDS(gene_list, "gene_listPrecursortoNeuroblast.rds")
gene_list <- readRDS("gene_listPrecursortoNeuroblast.rds")

gene_list <- row.names(gene_list)
df <- data.frame(time = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])
df <- subset(df, rownames(df) %in% colnames(cds)) 
colnames(df) <-'time'
for ( i in 1:length(gene_list)){
  
  DE <- gene_list[i]
  ite <- cds[row.names(cds) %in% DE,]
  df1 <- data.frame('Expression' = c(as.vector(log(exprs(ite)+1))),
                    'Pseudotime'= df$time,
                    'Cell_ID' = rownames(df))
  term <- predict(loess(Expression ~ Pseudotime, data=df1, span=0.75)) 
  df <- cbind(df,term)
}
colnames(df)[2:(length(gene_list)+1)] <- gene_list
df <- df %>% 
  arrange(time)
saveRDS(df,"PrecursortoNeuroblast_time_exp.rds")

## convert the data frame to time_series
t <- t(df)
colnames(t) <- df$time
t <- t %>% 
  as.data.frame()
t$gene <- row.names(t)
t <- t[-1,] %>% 
  relocate(gene) 
saveRDS(t,"PrecursortoNeuroblast.time_series.rds")
write.csv(t,"PrecursortoNeuroblast.time_series.csv",row.names = FALSE) 

x <- apply(X = t[,-1], MARGIN = 1, FUN = max)
y <- apply(X = t[,-1], MARGIN = 1, FUN = min)
m <- apply(X = t[,-1], MARGIN = 1, FUN = mean)
v <- apply(X = t[,-1], MARGIN = 1, FUN = sd)
df.t <- data.frame(gene = t[,1],max=x,min=y, mean=m, sd=v)
df.t$range <- df.t$max-df.t$min
saveRDS(df.t,"PrecursortoNeuroblast.df.t.rds")

SD <- df.t$sd
hist(SD, breaks=60)
MEAN <- df.t$mean
hist(MEAN, breaks=60)
Range <- df.t$range
hist(Range, breaks=60)
Mean_SD <- df.t[, c("mean", "sd")]
Range_SD <- df.t[, c("range", "sd")]
Mean_Range <- df.t[, c("mean", "range")]
plot(Mean_Range)

##Test the cutoff by range
remain <- subset(df.t, df.t$range>0.3) #2,734 genes
remain <- remain %>% arrange(range)

#Create time series on the "remain" genes
df <- data.frame(time = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])
df <- df[rownames(df)%in% colnames(cds),] %>% 
  as.data.frame()
colnames(df) <-'time'
len<- (max(df$time)-min(df$time))/(length(unique(df$time))-1)
new_time <- data.frame(time=seq(min(df$time),max(df$time),by=len))

gene_list <- remain$gene

for (i in 1:length(gene_list)){
  DE <- gene_list[i]
  ite <- cds[row.names(cds) %in% DE,]
  df1 <- data.frame('Expression' = c(as.vector(log(exprs(ite)+1))),
                    'Pseudotime'= df$time)
  term <- predict(loess(Expression ~ Pseudotime, data=df1, span=0.75),data.frame(Pseudotime=seq(min(df$time),max(df$time),by=len)))
  new_time <- cbind(new_time,term)
}
colnames(new_time)[2:ncol(new_time)] <- gene_list
new_time <- new_time %>% 
  arrange(time)
saveRDS(new_time,"PrecursortoNeuroblast_final_time_exp.rds")

tt <- scale(new_time[,-1])
new_time[,2:ncol(new_time)] <- tt
## it's log(exprs+1) then normalized
saveRDS(new_time,"PrecursortoNeuroblast_final_time_expN.rds")
write.csv(new_time,"PrecursortoNeuroblast_final_time_expN.csv",row.names = FALSE)






