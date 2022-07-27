library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ecodist)

seed_data <- read_csv("./data/drought_seed_production_data.csv")
nls_boot_pairs<- read.csv("./output/nls_boot_pairs_1000.csv")
pairs_abc <- c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
               "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
               "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
               "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
               "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
               "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


#get all the pairs where the fd is above 1 in the wet plots----
nls_boot_pairs$fd_superior <- ifelse(nls_boot_pairs$fd < 1, 1/nls_boot_pairs$fd, nls_boot_pairs$fd)
nls_boot_pairs$fd_sup_sp <- ifelse(nls_boot_pairs$fd <= 1, 1, 2)
nls_boot_pairs_sup <- nls_boot_pairs %>%
  filter(fd_sup_sp == 2 )
W_superior <- with(nls_boot_pairs_sup, sp_pair[treat=="W"])

boots_pairs_w_sup <- nls_boot_pairs %>%
  filter(sp_pair %in% W_superior)
boots_pairs_w_sup$treat <- factor(boots_pairs_w_sup$treat, levels = c("W", "D"))

nls_boot_pairs_unique <-  nls_boot_pairs %>%
  filter(sp_pair %in% pairs_abc)
nls_boot_pairs_unique$treat <- factor(nls_boot_pairs_unique$treat, levels = c("W", "D"))

nls_boot_pairs_unique$fd_superior <- ifelse(nls_boot_pairs_unique$fd < 1, 1/nls_boot_pairs_unique$fd, nls_boot_pairs_unique$fd)

nls_boot_pairs_unique$coexist <- ifelse((nls_boot_pairs_unique$snd > (1-1/nls_boot_pairs_unique$fd_superior)), 1, 0 )


boots_pairs_w_sup$treat <- factor(boots_pairs_w_sup$treat, levels = c("W", "D"))
boots_pairs_w_sup$label <- paste0(substr(boots_pairs_w_sup$focal, 1, 2), "-", substr(boots_pairs_w_sup$competitor, 1, 2))


# min/max fitness difference that permits coexistence
niche_differentiation <- seq(from = -.25, to = 1, by = 0.001)
niche_overlap <- 1-niche_differentiation
fitness_ratio_min <- niche_overlap
fitness_ratio_max <- 1/niche_overlap

df <- data.frame(niche_diff = niche_differentiation,
                 min_fitness_ratio = fitness_ratio_min,
                 max_fitness_ratio = fitness_ratio_max)
head(df)

#Text size issue
geom.text.size = 6
label.size = 3
element.size = (14/5)*label.size
theme.size = (14/5) * geom.text.size

#Coexistence plot
ggplot(boots_pairs_w_sup) + 
  theme_classic(base_size = 20)+
  coord_cartesian(xlim = c(-.05, 1), ylim = c(0.25, 10)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  annotate( "text", x = 0.8, y = .71, label = "Coexistence: \u03c1 < (Ki/Kj)", size = geom.text.size, fontface = 2) +
  annotate("text", x = 0.13, y = 6.5, label = "Competitive Exclusion: \u03c1 > (Ki/Kj)", size = geom.text.size, fontface = 2) +
  annotate("text", x = 0.13, y = .55, label = "Competitive Exclusion: \u03c1 > (Ki/Kj)", size = geom.text.size, fontface = 2) +
  geom_point(aes(x = snd, y = fd, color = treat, shape = treat, group = sp_pair), size = 4, stroke = 1) +
  scale_color_manual(values=c("W" ="#4E84C4", "D" = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c("W" = 15, "D" = 19), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  geom_text_repel(data = boots_pairs_w_sup %>%
                    filter(treat == "D"), 
                  aes(x = snd, y = fd, label = label, size = label.size), 
                  nudge_x = c(0.01, .025, -.025, 0.025, .025, 0,.025, -0.025, -.025, -0.025, -.025,-0.025,0.025,0.025,0.025 ), 
                  nudge_y = c(0.01, 0, 0, 0, 0,0.05, 0, 0, 0, 0, -0.01, 0,0,0,0 )) +
  geom_path(aes(x = snd, y = fd, group = sp_pair), arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(expand= c(.1,0)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)", size = 18) +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
          legend.title = element_text(size = theme.size),
          legend.text = element_text(size = theme.size),
          legend.position = c(.2, .025),
          legend.justification = c("bottom"),
          legend.box.just = "bottom",
          legend.margin = margin(1, 1, 1, 1),
          legend.direction = "horizontal") +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , size = F)+
  NULL
ggsave("./figures/boots_w_sup.jpeg", width = 10, height = 6)

#lambda
seed_data$treat <- factor(seed_data$treat, levels = c("W", "D"))

seed_data %>%
  filter(num_comp == 0) %>%
  ggplot(aes(x = focal, y = num_seeds, fill = factor(treat))) +
  theme_classic(base_size = 20) +
  theme(text = element_text(size = 18), 
    #axis.title = element_text(size = 18),
         #legend.title = element_text(size = 18),
         #legend.text = element_text(size = 15), 
         legend.position = c(.185, .11),
         legend.justification = c("bottom"),
         legend.box.just = "bottom",
         legend.margin = margin(1, 1, 1, 1)) +
  guides(fill = guide_legend(title = "Treatment:", title.position = "left"), col = guide_legend(nrow = 1), size = F )+
  geom_boxplot(position = "dodge") +
  scale_fill_manual(values=c("W" ="#4E84C4", "D" = "#D16103"), name = "treatment", labels = c("Ambient", "Reduced Rain")) +
  geom_text(aes(x = 1, y = 10000, label = "ns"), size = geom.text.size) +
  geom_text(aes(x = 2, y = 5000, label = "***"), size = geom.text.size) +
  geom_text(aes(x = 3, y = 3000, label = "ns"), size = geom.text.size) +
  geom_text(aes(x = 4, y = 3000, label = "ns"), size = geom.text.size) +
  geom_text(aes(x = 5, y = 4000, label = "ns"), size = geom.text.size) +
  geom_text(aes(x = 6, y = 1400, label = "***"), size = geom.text.size) +
  annotate("text", x = 1.8, y = 65, label = paste("*** = p < 0.001, ns = Not significant"), size = geom.text.size) +
  xlab(" ")+
  ylab("Fecundity (log scale)")+
  labs( color = "treatment") +
  scale_y_log10() +
  NULL
ggsave("./figures/real_lambda_log.pdf", width = 10, height = 6)


#TRAITS figure
  

  df_compile<- read_csv("./data/pca_mantel.csv")

  #these are from pca_traits.R, they provide r2 and p values
  ecodist :: mantel(dist_fd_d ~ dist_pc1)
  ecodist :: mantel(dist_fd_w ~ dist_pc1)
  ecodist :: mantel(dist_snd_d ~ dist_pc1)
  ecodist :: mantel(dist_snd_w ~ dist_pc1)
  ecodist :: mantel(dist_fd_diffs_abs ~ dist_pc1)
  ecodist :: mantel(dist_fd_diffs ~ dist_pc1)
  ecodist :: mantel(dist_snd_diffs ~ dist_pc1)
  
 fd_plot <- ggplot(df_compile)+
    geom_point(aes(x = pc1, y = fd_w ), color = "#4E84C4") +
    geom_point(aes(x = pc1, y = fd_d), color = "#D16103") +
    theme_classic() +
    annotate("text", x = 0.5, y = 6.6, label =paste("R^2 == ", "0.79"), size = 4, color = "#D16103", parse = T, hjust = 0) +
    annotate("text", x = 0.5, y = 6, label = "p = 0.003**", size = 4, color = "#D16103", hjust = 0) +
    annotate("text", x = 4, y = 1.3, label = paste("R^2 == ","0.44") , size = 4, color = "#4E84C4", parse = T, hjust = 0) +
    annotate("text", x = 4, y = 0.75, label = "p = .09", size = 4, color = "#4E84C4", hjust = 0) +
    #annotate("text", x= 4.8, y = 6.75, label = "A", size = 4, fontface = 2) + 
    stat_smooth(aes(x = pc1, y = fd_w), color = "#4E84C4", method = "lm", se=F) +
    stat_smooth(aes(x = pc1, y = fd_d), color = "#D16103", method = "lm", se = F) +
    ylab("Fitness differneces") +
    xlab("functional trait distances \n between pairs on pc1 axis")
  
  snd_plot <- ggplot(df_compile)+
    geom_point(aes(x = pc1, y = snd_w ), color = "#4E84C4") +
    geom_point(aes(x = pc1, y = snd_d), color = "#D16103") +
    theme_classic() +
    annotate("text", x = 0.5, y = 1.1, label =paste("R^2 == ", "0.075"), size = 4, color = "#D16103", parse = T, hjust = 0) +
    annotate("text", x = 0.5, y = 1, label = "p = 0.78", size = 4, color = "#D16103", hjust = 0) +
    annotate("text", x = 4, y = 0.14, label = paste("R^2 ==","-0.049"), size = 4, color = "#4E84C4", hjust = 0, parse = T) +  
    annotate("text", x = 4, y = 0.02, label = "p = 0.84", size = 4, color = "#4E84C4", hjust = 0) +
    #annotate("text", x= 4.75, y = 1, label = "B", size = 4, fontface = 2)+
    stat_smooth(aes(x = pc1, y = snd_w), color = "#4E84C4", method = "lm", se = F) +
    stat_smooth(aes(x = pc1, y = snd_d), color = "#D16103", method = "lm", se = F) +
    ylab("Stabilizing niche differneces") +
    xlab("functional trait distances \n between pairs on pc1 axis")
  
  snd_diff_plot <- ggplot(df_compile, aes(x = pc1, y = snd_diffs))+
    geom_point()+
    theme_classic() +
    geom_text_repel(aes(label = label), size = 3,
                    nudge_y = c(0.01, 0.01, -0.01, 0.01, -.01, .01 , .01, 0.01, .01, -.01, .01, .01, .01, .01, .01), 
                    nudge_x = c(.1, .1, .1, .1, .1, .1 , .1, 0.1, .1, .1, .1, .1, 0.1, 0.1, .1)) +
    annotate("text", x = 0.5, y = .68, label = paste("R^2 ==", "0.23"), hjust = 0, parse = T) +
    annotate("text", x = 0.5, y = .61, label = "p = 0.41", size = 4, hjust = 0) +
    #annotate("text", x= 4.8, y = 0.68, label = "D", size = 4, fontface = 2)+
    stat_smooth(method = "lm", color = "black", se = F) +
    ylab("Change in stabilizing niche \n differneces between treatments") +
    xlab("functional trait distances \n between pairs on pc1 axis")
  
  fd_diff_plot_abs <- ggplot(df_compile, aes(x = pc1, y = fd_diff_abs))+
    geom_point()+
    theme_classic() +
    geom_text_repel(aes(label = label), size = 3,
                    nudge_y = c(.1, -.1, -.01, .1, .1, .1 , .1, 0, .1, -.1, .1, .1, .1, .05, .1), 
                    nudge_x = c(.1, .1, .1, .1, .1, .1 , .1, 0, .1, .1, .1, .1, 0.1, 0, .1))+
    annotate("text", x = 0.5, y = 4, label = paste("R^2 ==","0.43"), size = 4, parse = T, hjust = 0)+
    annotate("text",x = 0.5, y = 3.6, label = "p = 0.021*", size = 4, hjust = 0) +
    #annotate("text", x= 4.75, y = 4, label = "C", size = 4, fontface = 2)+
    stat_smooth(method = "lm", color = "black", se= F) +
    ylab("Change in fitness differneces \n between treatments (absolute value)") +
    xlab("functional trait distances \n between pairs on pc1 axis")
  
diff_plots <- ggarrange(fd_diff_plot_abs, snd_diff_plot, labels = c("C", "D"))
fd_snd_plot <- ggarrange(fd_plot, snd_plot, labels = c("A", "B"))
ggarrange( fd_snd_plot, diff_plots, nrow = 2)
  
ggsave("./figures/trait_diff.jpeg", width = 10, height = 6)

#Supplemental Coexistence plot

ggplot(boots_pairs_w_sup) + 
  theme_classic(base_size = 20)+
  coord_cartesian(xlim = c(-.05, 1), ylim = c(0.1, 15)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treat, shape = treat, group = sp_pair), size = 2, stroke = 1) +
  geom_errorbar(aes(x = snd, ymin = (fd_low), ymax = (fd_high))) +
  geom_errorbarh(aes(y = fd, xmin = (snd_low), xmax = (snd_high))) +
  scale_color_manual(values=c("W" ="#4E84C4", "D" = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c("W" = 15, "D" = 16), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  #geom_path(aes(x = snd, y = fd, group = sp_pair), arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  #scale_x_continuous(expand= c(.1,0)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)", size = 18) +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
        legend.title = element_text(size = theme.size),
        legend.text = element_text(size = theme.size),
        legend.position = c(0.88, .035),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical", 
        strip.background = element_blank()) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , size = F)+
 facet_wrap(vars(label), nrow = 4, scales ="free")+
  NULL
ggsave("./figures/coexistence_facet.jpeg", width = 12, height = 10)
  