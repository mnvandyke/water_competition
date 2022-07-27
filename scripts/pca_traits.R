#Trait PCA 2019
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggfortify)
library(vegan)
library(ggpubr)
library(ecodist)

gkt <- read.csv("./data/traits_gk.csv")
gkt[23, "seed_mass_g"] <- 0.01396391 #seed weight for HOMU missing - this is from purity data (includes wings)
gkt<-gkt[-21, ] # BRMA - no data
gkt <- dplyr::select(gkt, !c(leaf_pH, lar_cm2_g, seed_size_mm3, lai_la_canopy, d15N, CN_ratio))
rownames(gkt) <- gkt$species

nls_boot_pairs <- read.csv("./output/nls_boot_pairs_1000.csv")
nls_boot_pairs$fd_superior <- ifelse(nls_boot_pairs$fd < 1, 1/nls_boot_pairs$fd, nls_boot_pairs$fd)

gkt_pca <- prcomp(gkt[,2:12], scale = T)
summary(gkt_pca)
autoplot(gkt_pca, data = gkt, 
         loadings = TRUE, 
         loadings.colour = 'red', 
         loadings.label = T, 
         loadings.label.size = 3, 
         loadings.label.vjust = 1.2, label = T)
ggsave("./figures/pca.pdf")

gkt_pca$x
gkt_pca$x[18,1]
screeplot(gkt_pca)

pc1<- c(gkt_pca$x[15,1], gkt_pca$x[24,1], gkt_pca$x[22,1], gkt_pca$x[18,1], gkt_pca$x[19,1], gkt_pca$x[1,1])
pc2<- c(gkt_pca$x[15,2], gkt_pca$x[24,2], gkt_pca$x[22,2], gkt_pca$x[18,2], gkt_pca$x[19,2], gkt_pca$x[1,2])
#Dist matrices for fd and snd
fd_d  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd <- nls_boot_pairs %>% 
      filter(treat == "D" & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_superior) %>%
      unlist()
    fd_d[sp1, sp2] <- current_fd  }}
dist_fd_d<- as.dist(fd_d)

snd_d <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(snd_d)) {  
  for(sp2 in colnames(snd_d)) {    
    current_snd <- nls_boot_pairs %>% 
      filter(treat == "D" & focal == sp1 & competitor == sp2) %>% 
      select(snd) %>%
      unlist()
    snd_d[sp1, sp2] <- current_snd  }}

dist_snd_d <- as.dist(snd_d)

fd_w  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_w)) {  
  for(sp2 in colnames(fd_w)) {    
    current_fd <- nls_boot_pairs %>% 
      filter(treat == "W" & focal == sp1 & competitor == sp2) %>% 
      select(fd_superior) %>%
      unlist()
    fd_w[sp1, sp2] <- current_fd  }}
dist_fd_w<- as.dist(fd_w)

snd_w  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(snd_w)) {  
  for(sp2 in colnames(snd_w)) {    
    current_snd <- nls_boot_pairs %>% 
      filter(treat == "W" & focal == sp1 & competitor == sp2) %>% 
      select(snd) %>%
      unlist()
    snd_w[sp1, sp2] <- current_snd  }}
dist_snd_w <- as.dist(snd_w)

dist_pc1<- dist(pc1)
dist_pc2 <- dist(pc2)
vegan :: mantel(dist_pc1, dist_snd_w)
vegan ::mantel(dist_pc1, dist_fd_w)
vegan ::mantel(dist_pc1, dist_snd_d)
vegan :: mantel(dist_pc1, dist_fd_d)
ecodist :: mantel(dist_fd_d ~ dist_pc1)
ecodist :: mantel(dist_fd_w ~ dist_pc1)
ecodist :: mantel(dist_snd_w ~ dist_pc1)
ecodist :: mantel(dist_snd_d ~ dist_pc1)

pc1_vec <- as.vector(dist_pc1)
fd_d_vec <- as.vector(dist_fd_d)
plot(pc1_vec, fd_d_vec)
abline(1,1)

vegan ::mantel(dist_pc2, dist_snd_w)
vegan ::mantel(dist_pc2, dist_fd_w)
vegan ::mantel(dist_pc2, dist_snd_d)
vegan ::mantel(dist_pc2, dist_fd_d)

#Do traits dissimilarities between pairs correlate with magnitude that their snd and fd change between treatments
diffs <- nls_boot_pairs
diffs$fd_diff_abs = abs(diffs$fd[diffs$treat == "W"] - diffs$fd[diffs$treat == "D"])
diffs$snd_diff = abs(diffs$snd[diffs$treat == "W"] - diffs$snd[diffs$treat == "D"])
diffs$fd_diff = diffs$fd[diffs$treat == "W"] - diffs$fd[diffs$treat == "D"]

diffs <- diffs %>% 
  dplyr :: select(sp_pair, focal, competitor, treat,fd_diff, fd_diff_abs, snd_diff) %>%
  filter(treat == "W") 
diffs$treat <- NULL

# snd and fd diff distance matrices
#make a matrix of snd differences between treatments for each pair
snd_diffs <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(snd_diffs)) {  
  for(sp2 in colnames(snd_diffs)) {    
    snd_change <- diffs %>% 
      filter(focal == sp1 & competitor == sp2) %>% 
      dplyr::select(snd_diff) %>%
      unlist()
    snd_diffs[sp1, sp2] <- snd_change  }}
dist_snd_diffs<- as.dist(snd_diffs)
snd_distances_vec <- as.vector(dist_snd_diffs)


fd_diffs  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_diffs)) {  
  for(sp2 in colnames(fd_diffs)) {    
    fd_change <- diffs %>% 
      filter(focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_diff) %>%
      unlist()
    fd_diffs[sp1, sp2] <- fd_change  }}
dist_fd_diffs<- as.dist(fd_diffs)
fd_diffs

fd_diffs_abs  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_diffs_abs)) {  
  for(sp2 in colnames(fd_diffs_abs)) {    
    fd_change <- diffs %>% 
      filter(focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_diff_abs) %>%
      unlist()
    fd_diffs_abs[sp1, sp2] <- fd_change  }}
dist_fd_diffs_abs<- dist(fd_diffs_abs)


vegan ::mantel(dist_pc1, dist_snd_diffs)
vegan :: mantel(dist_pc1, dist_fd_diffs)
vegan ::mantel(dist_pc1, dist_fd_diffs_abs) # pca trait differences correlated with the amount fitness diffs change
#make them into vectors and plot- check over this whole thing
pc1_vec <- as.vector(dist_pc1)
fd_diffs_abs_vec <- as.vector(dist_fd_diffs_abs)
fd_diffs_vec <- as.vector(dist_fd_diffs)

plot(pc1_vec, fd_diffs_vec)
abline(1,-1)

vegan ::mantel(dist_pc2, dist_snd_diffs)
vegan ::mantel(dist_pc2, dist_fd_diffs)
vegan ::mantel(dist_pc2, dist_fd_diffs_abs)

vegan :: mantel(dist_pc1, dist_fd_diffs_abs)
ecodist :: mantel(dist_fd_diffs_abs ~ dist_pc1)
ecodist :: mantel(dist_fd_diffs ~ dist_pc1)

ggplot( x = pc1_vec, y = fd_diffs_abs_vec) +
  geom_point(aes())

#Combine all to a data frame

fd_d_df <- as.data.frame(gdata :: unmatrix(fd_d)) 
colnames(fd_d_df) = "fd_d"
fd_d_df <- tibble::rownames_to_column(fd_d_df, "sp_pair")

fd_w_df <- as.data.frame(gdata :: unmatrix(fd_w))
colnames(fd_w_df) = "fd_w"
fd_w_df <- tibble::rownames_to_column(fd_w_df, "sp_pair")

snd_w_df <- as.data.frame(gdata :: unmatrix(snd_w))
colnames(snd_w_df) = "snd_w"
snd_w_df <- tibble::rownames_to_column(snd_w_df, "sp_pair")


snd_d_df <- as.data.frame(gdata :: unmatrix(snd_d))
colnames(snd_d_df) = "snd_d"
snd_d_df<- tibble::rownames_to_column(snd_d_df, "sp_pair")

fd_diff_df <-as.data.frame(gdata :: unmatrix(fd_diffs))
colnames(fd_diff_df) = "fd_diff"
fd_diff_df <- tibble::rownames_to_column(fd_diff_df, "sp_pair")

fd_diffs_abs_df <- as.data.frame(gdata :: unmatrix(fd_diffs_abs))
colnames(fd_diffs_abs_df) = "fd_diff_abs"
fd_diffs_abs_df <- tibble::rownames_to_column(fd_diffs_abs_df, "sp_pair")

snd_diffs_df <- as.data.frame(gdata :: unmatrix(snd_diffs))
colnames(snd_diffs_df) = "snd_diffs"
snd_diffs_df<- tibble::rownames_to_column(snd_diffs_df, "sp_pair")

pc1_mat <- as.matrix(dist(pc1, upper = T))
rownames(pc1_mat)<- c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
colnames(pc1_mat)<- c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
pc1_df <- as.data.frame(gdata::unmatrix(pc1_mat))
colnames(pc1_df) = "pc1"
pc1_df<- tibble::rownames_to_column(pc1_df, "sp_pair")


df_compile <- left_join(fd_d_df, fd_w_df, by = "sp_pair") %>%
  left_join(.,snd_d_df, by = "sp_pair") %>%
  left_join(.,snd_w_df, by = "sp_pair") %>%
  left_join(.,fd_diff_df, by = "sp_pair") %>%
  left_join(.,fd_diffs_abs_df, by = "sp_pair") %>%
  left_join(.,snd_diffs_df, by = "sp_pair") %>%
  left_join(., pc1_df, by = "sp_pair")

df_compile <- df_compile[-c(1, 8, 15, 22, 29, 36), ]

w_pairs<- c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", "PLER:FEMI", "SACO:FEMI", "URLI:FEMI",
"ACWR:HOMU", "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile<- df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label <- paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel.csv")

ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) +
  geom_point() +
  geom_smooth(method = "lm")


## FIGURES##

fd_d_plot <-
  ggplot(df_compile, aes(x = pc1, y = fd_d))+
  geom_point()+
  geom_text_repel(aes(label = label), size = 3,
                  nudge_y = c(.1, .1, -.01, .1, .1, .1 , .1, .1, .1, .1, .1, 0, .1, .1, .1), 
                  nudge_x = c(.1, -.1, .1, .1, .1, .1 , .1, .1, .1, .1, .1, .1, .1, .1, .1))+
  annotate("text", x = 1.2, y = 6.5, label =paste("R^2 == ", "0.79"), size = 4, parse = T) +
  annotate("text", x = 1.2, y = 6, label = "p = 0.004", size = 4) +
  stat_smooth(method = "lm")+
  theme_classic() +
  ylab("Fitness differneces in \n reduced rain plots") +
  xlab("functional trait distances \n between pairs on pc1 axis")


fd_w_plot <- 
  ggplot(df_compile, aes(x = pc1, y = fd_w))+
  geom_point() +
  theme_classic() +
  geom_text_repel(aes(label = label), size = 3,
                  nudge_y = c(.1, .1, -.01, .1, .1, .1 , .1, .1, .1, .1, .1, 0, .1, .1, .1), 
                  nudge_x = c(.1, -.1, .1, .1, .1, .1 , .1, .1, .1, .1, .1, .1, .1, .1, .1))+
  annotate("text", x = 1, y = 5, label = "R = 0.44\np = .08", size = 4, fontface = 2) +
  stat_smooth(method = "lm") +
  ylab("Fitness differneces in \n ambient plots") +
  xlab("functional trait distances \n between pairs on pc1 axis")


fd_plots <- ggarrange(fd_w_plot, fd_d_plot)
fd_diff_plot <- ggplot(df_compile, aes(x = pc1, y = fd_diff))+
  geom_point()+
  theme_classic() +
  geom_text_repel(aes(label = label), size = 3,
                  nudge_y = c(.1, .1, -.01, .1, .1, -.1 , -.1, .1, .1, .1, .1, 0, .1, .1, .1), 
                  nudge_x = c(.1, -.1, .1, .1, .1, 0, .1, .1, .1, .1, .1, .1, .1, .1, .1))+
  annotate("text", x = 4.5, y = 2.5, label = as.character(expression(paste(R^{2}*" = -0.66\np = .005"))), size = 4, fontface = 2, parse = T) +
  stat_smooth(method = "lm") +
  ylab("Change in fitness differneces \n between treatments") +
  xlab("functional trait distances \n between pairs on pc1 axis")

fd_diff_plot_abs <- ggplot(df_compile, aes(x = pc1, y = fd_diff_abs))+
  geom_point()+
  theme_classic() +
  geom_text_repel(aes(label = label), size = 3,
                  nudge_y = c(.1, -.1, -.01, .1, .1, .1 , .1, 0, .1, -.1, .1, .1, .1, .05, .1), 
                  nudge_x = c(.1, .1, .1, .1, .1, .1 , .1, 0, .1, .1, .1, .1, 0, 0, .1))+
  annotate("text", x = 1, y = 3.5, label = "R = 0.43\np = .017 ***", size = 4, fontface = 2) +
  stat_smooth(method = "lm") +
  ylab("Change in fitness differneces between \n treatments (absolute value)") +
  xlab("functional trait distances \n between pairs on pc1 axis")

fd_diff_plots <- ggarrange(fd_diff_plot_abs,fd_diff_plot)

ggarrange(fd_plots, fd_diff_plots,  nrow = 2)  

ggarrange(fd_diff_plot_abs, fd_diff_plot, nrow = 2)


fd_plot <- ggplot(df_compile)+
  geom_point(aes(x = pc1, y = fd_w ), color = "dodgerblue2") +
  geom_point(aes(x = pc1, y = fd_d), color = "#D55E00") +
  theme_classic() +
  annotate("text", x = 1.2, y = 6.5, label =paste("R^2 == ", "0.79"), size = 4, color = "#D55E00", parse = T) +
  annotate("text", x = 1.2, y = 6, label = "p = 0.004***", size = 4, color = "#D55E00") +
  annotate("text", x = 4.5, y = 1, label = "R = 0.44\np = .08", size = 4, color = "dodgerblue2") +
  stat_smooth(aes(x = pc1, y = fd_w), color = "dodgerblue2", method = "lm") +
  stat_smooth(aes(x = pc1, y = fd_d), color = "#D55E00", method = "lm") +
  ylab("Fitness differneces") +
  xlab("functional trait distances \n between pairs on pc1 axis")

snd_plot <- ggplot(df_compile)+
  geom_point(aes(x = pc1, y = snd_w ), color = "dodgerblue2") +
  geom_point(aes(x = pc1, y = snd_d), color = "#D55E00") +
  theme_classic() +
  annotate("text", x = 1.2, y = 1.2, label =paste("R^2 == ", "-0.049"), size = 4, color = "#D55E00", parse = T) +
  annotate("text", x = 1.2, y = 1, label = "p = 0.84", size = 4, color = "#D55E00") +
  annotate("text", x = 4.5, y = 0.125, label = "R = 0.075\np = .78", size = 4, color = "dodgerblue2") +  stat_smooth(aes(x = pc1, y = snd_w), color = "dodgerblue2", method = "lm") +
  stat_smooth(aes(x = pc1, y = snd_d), color = "#D55E00", method = "lm") +
  ylab("Stabilizing niche differneces") +
  xlab("functional trait distances \n between pairs on pc1 axis")

snd_diff_plot <- ggplot(df_compile, aes(x = pc1, y = snd_diffs))+
  geom_point()+
  theme_classic() +
  geom_text_repel(aes(label = label), size = 3,
                  nudge_y = c(0.01, 0.01, 0.01, 0.01, -.01, .01 , .01, 0.01, .01, -.01, .01, .01, .01, .01, .01), 
                  nudge_x = c(.1, .1, .1, .1, .1, .1 , .1, 0.1, .1, .1, .1, .1, 0.1, 0.1, .1))+
  annotate("text", x = 1, y = .6, label = "R = 0.23\np = 0.41", size = 4, fontface = 2) +
  stat_smooth(method = "lm") +
  ylab("Change in stabilizing niche \n differneces between treatments") +
  xlab("functional trait distances \n between pairs on pc1 axis")

diff_plots <- ggarrange(fd_diff_plot_abs, snd_diff_plot)
fd_snd_plot <- ggarrange(fd_plot, snd_plot)
ggarrange( fd_snd_plot, diff_plots, nrow = 2)

