</> R
## R version 4.1.2 (2021-11-01)

library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggeasy)
library(gridExtra)

#Hierachy clustering
new_time <- read.csv("PrecursortoNeuroblast_final_time_expN.csv", check.names = FALSE)
time <- new_time[,1]
gene_ts <- t(new_time[,-1]) %>% 
  as.data.frame()
hclust_res <- hclust(dist(gene_ts), method = "ward.D2")
ss_cluster <- cutree(hclust_res, k = 13)
ss_cluster <- as.data.frame(ss_cluster)
ss_cluster$gene <- rownames(ss_cluster)

saveRDS(ss_cluster, "ss_cluster_PrecursortoNeuroblast.rds")

ts <- read.csv("PrecursortoNeuroblast_final_time_expN.csv",check.names = FALSE) %>%
  t() %>%
  as.data.frame()
colnames(ts) <- ts[1,]
mean_df <- data.frame(time=as.numeric(colnames(ts)))
ts <- ts[-1,]
ts$gene <- rownames(ts)

result <- ss_cluster
ts_normal <- left_join(result,ts)

for (i in 1:length(unique(ts_normal$ss_cluster))){
  write.csv(ts_normal[ts_normal$ss_cluster==i,],paste("Cluster",i,".csv",sep = ""),row.names = F)
}

p <- list()
for (i in 1:length(unique(ts_normal$ss_cluster))){
  
  
  clr <- read.csv(paste("Cluster",i,".csv",sep = ""),check.names = F)
  clr_ts <- clr[,-c(1,2)] %>%
    sapply(as.numeric) %>%
    colMeans()
  mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
  mean_df <- cbind(mean_df,mean_ts$mean)
  colnames(mean_df)[i+1] <- paste0('cluster',i,sep="")
  #clr <- rbind(clr,c("Redline","mean",clr_ts))
  clr_gg<-clr %>%
    gather(time, exprs,-ss_cluster,-gene)
  
  clr_gg$time <- as.numeric(clr_gg$time)
  clr_gg$exprs <- as.numeric(clr_gg$exprs)
  
  #cols <- c("gray","red")
  
  
  p[[i]]<-ggplot() +
    geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
    geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
    xlab('pseudotime') +
    ylab('log(exprs+1)')+
    #ylim(0,5)+
    ggtitle(paste("Cluster",i,sep = ""))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggeasy::easy_center_title()
  
}

mean_df <- t(mean_df) %>%
  as.data.frame()

colnames(mean_df) <- mean_df[1,]
mean_df <- mean_df[-1,]
mean_df$cluster <- rownames(mean_df)

mean_gg<-mean_df %>%
  gather(time, exprs,-cluster)

mean_gg$time <- as.numeric(mean_gg$time)
mean_gg$exprs <- as.numeric(mean_gg$exprs)

pdf("cluster_PrecursortoNeuroblast.pdf", width = 12, height = 4)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]], p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],
             p[[11]], p[[12]],p[[13]],nrow = 4)
dev.off()

## Combine clusters after reading GO and inspect patterns

ts <- read.csv("PrecursortoNeuroblast_final_time_expN.csv",check.names = FALSE) %>%
  t() %>%
  as.data.frame()

colnames(ts) <- ts[1,]
mean_df <- data.frame(time=as.numeric(colnames(ts)))
ts <- ts[-1,]
ts$gene <- rownames(ts)

result <- ss_cluster
ts_normal <- left_join(result,ts)

write.csv(ts_normal[ts_normal$ss_cluster %in% c(1,3,5,11),],paste("Cluster",135.11,".csv",sep = ""),row.names = FALSE)
write.csv(ts_normal[ts_normal$ss_cluster %in% c(4,7),],paste("Cluster",47,".csv",sep = ""),row.names = FALSE)
write.csv(ts_normal[ts_normal$ss_cluster %in% c(2,6,8),],paste("Cluster",268,".csv",sep = ""),row.names = FALSE)
write.csv(ts_normal[ts_normal$ss_cluster %in% c(10,12),],paste("Cluster",10.12,".csv",sep = ""),row.names = FALSE)

saveRDS(ts_normal[ts_normal$ss_cluster %in% c(1,3,5,11),],paste("Cluster",135.11,".rds"))
saveRDS(ts_normal[ts_normal$ss_cluster %in% c(4,7),],paste("Cluster",47,".rds"))
saveRDS(ts_normal[ts_normal$ss_cluster %in% c(2,6,8),],paste("Cluster",268,".rds"))
saveRDS(ts_normal[ts_normal$ss_cluster %in% c(10,12),],paste("Cluster",10.12,".rds"))
saveRDS(ts_normal[ts_normal$ss_cluster %in% c(9),],paste("Cluster",9,".rds"))
saveRDS(ts_normal[ts_normal$ss_cluster %in% c(13),],paste("Cluster",13,".rds"))

p <- list()
clr <- read.csv(paste("Cluster",135.11,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[2] <- paste0('cluster',135.11,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[1]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",135.11,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

clr <- read.csv(paste("Cluster",47,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[3] <- paste0('cluster',47,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[2]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",47,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

clr <- read.csv(paste("Cluster",268,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[4] <- paste0('cluster',268,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[3]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",268,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

clr <- read.csv(paste("Cluster",10.12,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[5] <- paste0('cluster',10.12,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[4]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",10.12,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

clr <- read.csv(paste("Cluster",9,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[6] <- paste0('cluster',9,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[5]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",9,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

clr <- read.csv(paste("Cluster",13,".csv",sep = ""),check.names = F)
clr_ts <- clr[,-c(1,2)] %>%
  sapply(as.numeric) %>%
  colMeans()
mean_ts <- data.frame(time=as.numeric(colnames(clr)[3:ncol(clr)]),mean=clr_ts)
mean_df <- cbind(mean_df,mean_ts$mean)
colnames(mean_df)[7] <- paste0('cluster',13,sep="")
clr_gg<-clr %>%
  gather(time, exprs,-ss_cluster,-gene)
clr_gg$time <- as.numeric(clr_gg$time)
clr_gg$exprs <- as.numeric(clr_gg$exprs)

p[[6]]<-ggplot() +
  geom_line(data=clr_gg,aes(x = time, y = exprs,group=gene),color="grey") +
  geom_line(data=mean_ts,aes(x=time,y = mean), color = "red",linewidth=1)+
  xlab('pseudotime') +
  ylab('log(exprs+1)')+
  #ylim(0,5)+
  ggtitle(paste("Cluster",13,sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()

mean_df <- t(mean_df) %>%
  as.data.frame()

colnames(mean_df) <- mean_df[1,]
mean_df <- mean_df[-1,]
mean_df$cluster <- rownames(mean_df)
mean_gg<-mean_df %>%
  gather(time, exprs,-cluster)
mean_gg$time <- as.numeric(mean_gg$time)
mean_gg$exprs <- as.numeric(mean_gg$exprs)

pdf("cluster_PrecursortoNeuroblast.combine.pdf", width = 12, height = 4)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],nrow = 1)
dev.off()

pdf("mean_PrecursortoNeuroblast.combine.pdf", width = 6, height = 4)
ggplot(mean_gg, aes(x = time, y = exprs, color = cluster)) +
  geom_line()+xlab('pseudotime') +
  ylab('log(exprs+1)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggeasy::easy_center_title()
dev.off()
