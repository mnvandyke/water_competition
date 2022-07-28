### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### this script creates a pca of functional traits from 23 species in the community (figure ED2)
### it then uses pca axes to compare trait differences with changes in stabilizing niche and fitness differences

library(dplyr)
library(tidyr)
library(tidyverse)
library(ggfortify)
library(vegan)
library(ggpubr)
library(gdata)


gkt <- read.csv("./data/traits_gk.csv")
rownames(gkt) <- gkt$species

nls_boot_pairs <- read.csv("./output/nls_boot_pairs_1000_full_model.csv") # or run model script first to get this
nls_boot_pairs$fd_superior <- ifelse(nls_boot_pairs$fd < 1, 1/nls_boot_pairs$fd, nls_boot_pairs$fd)

gkt_pca <- prcomp(gkt[,2:12], scale = T)
summary(gkt_pca)
autoplot(gkt_pca, data = gkt, 
         loadings = TRUE, 
         loadings.colour = 'red', 
         loadings.label = T, 
         loadings.label.size = 3, 
         loadings.label.vjust = 1.2, label = T)
#ggsave("./figures/pca.pdf")

gkt_pca$x
gkt_pca$x["ACWR",1]
screeplot(gkt_pca)

#Create PC vectors with the six species in the experiment
pc1<- c(gkt_pca$x["ACWR",1], gkt_pca$x["FEMI",1], gkt_pca$x["HOMU",1], gkt_pca$x["PLER",1], gkt_pca$x["SACO",1], gkt_pca$x["URLI",1])
pc2<- c(gkt_pca$x["ACWR",2], gkt_pca$x["FEMI",2], gkt_pca$x["HOMU",2], gkt_pca$x["PLER",2], gkt_pca$x["SACO",2], gkt_pca$x["URLI",2])
#Dist matrices for fd and snd
fd_d  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd <- nls_boot_pairs %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_superior) %>%
      unlist()
    fd_d[sp1, sp2] <- current_fd  }}
dist_fd_d<- as.dist(fd_d)

snd_d <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(snd_d)) {  
  for(sp2 in colnames(snd_d)) {    
    current_snd <- nls_boot_pairs %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      select(snd) %>%
      unlist()
    snd_d[sp1, sp2] <- current_snd  }}

dist_snd_d <- as.dist(snd_d)

fd_w  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(fd_w)) {  
  for(sp2 in colnames(fd_w)) {    
    current_fd <- nls_boot_pairs %>% 
      filter(treatment == 1 & focal == sp1 & competitor == sp2) %>% 
      select(fd_superior) %>%
      unlist()
    fd_w[sp1, sp2] <- current_fd  }}
dist_fd_w<- as.dist(fd_w)

snd_w  <- matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

for(sp1 in rownames(snd_w)) {  
  for(sp2 in colnames(snd_w)) {    
    current_snd <- nls_boot_pairs %>% 
      filter(treatment == 1 & focal == sp1 & competitor == sp2) %>% 
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
diffs$fd_diff_abs = abs(diffs$fd[diffs$treat == 1] - diffs$fd[diffs$treat == 2])
diffs$snd_diff = abs(diffs$snd[diffs$treat == 1] - diffs$snd[diffs$treat == 2])
diffs$fd_diff = diffs$fd[diffs$treat == 1] - diffs$fd[diffs$treat == 2]

diffs <- diffs %>% 
  dplyr :: select(sp_pair, focal, competitor, treatment,fd_diff, fd_diff_abs, snd_diff) %>%
  filter(treatment == 1) 


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
#make them into vectors and plot
pc1_vec <- as.vector(dist_pc1)
fd_diffs_abs_vec <- as.vector(dist_fd_diffs_abs)
fd_diffs_vec <- as.vector(dist_fd_diffs)

plot(pc1_vec, fd_diffs_vec)
abline(1,-1)

vegan ::mantel(dist_pc2, dist_snd_diffs)
vegan ::mantel(dist_pc2, dist_fd_diffs)
vegan ::mantel(dist_pc2, dist_fd_diffs_abs)
vegan :: mantel(dist_pc1, dist_fd_diffs_abs)

ggplot() +
  geom_point(aes(x = pc1_vec, y = fd_diffs_abs_vec))

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

df_compile <- df_compile %>%
  filter(fd_w != fd_d)

w_pairs<- c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", 
            "PLER:FEMI", "SACO:FEMI", "URLI:FEMI","ACWR:HOMU", 
            "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", 
            "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile<- df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label <- paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel.csv")

ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) +
  geom_point() +
  geom_smooth(method = "lm")

