library(ggplot2)
library(ggthemes)
# library(ggpubr)
ggplot2::theme_set(theme_bw())

res_path <- "/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_trajectory_fate_prob_moscot/prob_res"
plot_path <- "/Users/zhaojia/Documents/Yale/Research/Projects - Rui/20240304, ICNs - Yvonne/oin_trajectory_fate_prob_moscot/plots"


cutoff_1 <- 3.9659309738075756
cutoff_2 <- 12.815343374020156


### H-G
prob_organ_1 <- read.csv(paste0(res_path, "/prob_H_G_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_H_G_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("ENS_E14" = "#8299FF", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
pdf(paste0(plot_path, "/res_H_G.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("ENS_E14" = "#8299FF", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
pdf(paste0(plot_path, "/res_H_G_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



### H-L
prob_organ_1 <- read.csv(paste0(res_path, "/prob_H_L_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_H_L_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
pdf(paste0(plot_path, "/res_H_L.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
            aes(x = pseudotime1to4, 
                y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
pdf(paste0(plot_path, "/res_H_L_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



### H-P
prob_organ_1 <- read.csv(paste0(res_path, "/prob_H_P_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_H_P_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_pancrease" = "#80E380", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_H_P.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
            aes(x = pseudotime1to4, 
                y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_pancrease" = "#80E380", "ICN.E14.3" = "#FF9380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Heart") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_H_P_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



### G-L
prob_organ_1 <- read.csv(paste0(res_path, "/prob_G_L_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_G_L_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "ENS_E14" = "#8299FF")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Gut") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_G_L.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
            aes(x = pseudotime1to4, 
                y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "ENS_E14" = "#8299FF")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Gut") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_G_L_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



### G-P
prob_organ_1 <- read.csv(paste0(res_path, "/prob_G_P_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_G_P_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_pancrease" = "#80E380", "ENS_E14" = "#8299FF")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Gut") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_G_P.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
            aes(x = pseudotime1to4, 
                y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_pancrease" = "#80E380", "ENS_E14" = "#8299FF")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Gut") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_G_P_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



### L-P
prob_organ_1 <- read.csv(paste0(res_path, "/prob_L_P_organ1.csv"))
prob_organ_2 <- read.csv(paste0(res_path, "/prob_L_P_organ2.csv"))
prob <- rbind(prob_organ_1, prob_organ_2)
p <- ggplot(data = prob,
       aes(x = pseudotime1to4, 
           y = prob_organ1, color = orig.ident)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "E14.5_pancrease" = "#80E380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Lung") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_L_P.pdf"), width = 4.48, height = 3.23)
p
dev.off()

p <- ggplot(data = prob,
            aes(x = pseudotime1to4, 
                y = prob_organ1, color = orig.ident)) +
  # geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(aes(color = orig.ident), method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("E14.5_Lung" = "#DE9BF4", "E14.5_pancrease" = "#80E380")) +
  geom_vline(xintercept = cutoff_1, linewidth = 0.4) + 
  geom_vline(xintercept = cutoff_2, linewidth = 0.4) +
  # scale_color_viridis_d(option = "C") +
  #facet_grid(~ID) +
  theme(aspect.ratio = 1) +
  labs(color = "") +
  xlab("Pseudotime") +
  ylab("Fate to Lung") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
pdf(paste0(plot_path, "/res_L_P_lines.pdf"), width = 4.48, height = 3.23)
p
dev.off()



