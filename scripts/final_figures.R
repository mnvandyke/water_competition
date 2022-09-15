### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

##Reproduces most figures and tables in the paper, 
#uses output from most other scripts including:
#(nls_boots_all_model.R, pca_traits.R, lambda_plants.R, n_alpha_ratios.R, soil_gwc.R)
#structural analysis tables are in another script (structural_table.R)
#figure ED2 is created in pca_traits.R

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gt)
library(webshot)

options(scipen = 5)

seed_data <- read_csv("./data/drought_seed_production_data.csv")
seed_data$treatment <- ifelse(seed_data$treat == "W", 1, 2)

nls_boot_pairs<- read.csv("./output/nls_boot_pairs_1000_full_model.csv")

pairs_abc <- c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
               "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
               "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
               "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
               "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
               "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


#get all the pairs where the fitness difference is above 1 in the wet plots----
nls_boot_pairs$fd_superior <- ifelse(nls_boot_pairs$fd < 1, 1/nls_boot_pairs$fd, nls_boot_pairs$fd)
nls_boot_pairs$fd_sup_sp <- ifelse(nls_boot_pairs$fd <= 1, 1, 2)
nls_boot_pairs_sup <- nls_boot_pairs %>%
  filter(fd_sup_sp == 2 )
W_superior <- with(nls_boot_pairs_sup, sp_pair[treatment==1])

boots_pairs_w_sup <- nls_boot_pairs %>%
  filter(sp_pair %in% W_superior)
boots_pairs_w_sup$treat <- factor(boots_pairs_w_sup$treat, levels = c(1, 2))

nls_boot_pairs_unique <-  nls_boot_pairs %>%
  filter(sp_pair %in% pairs_abc)
nls_boot_pairs_unique$treatment <- factor(nls_boot_pairs_unique$treatment, levels = c(1, 2))

nls_boot_pairs_unique$fd_superior <- ifelse(nls_boot_pairs_unique$fd < 1, 1/nls_boot_pairs_unique$fd, nls_boot_pairs_unique$fd)

nls_boot_pairs_unique$coexist <- ifelse((nls_boot_pairs_unique$snd > (1-1/nls_boot_pairs_unique$fd_superior)), 1, 0 )


boots_pairs_w_sup$treatment <- factor(boots_pairs_w_sup$treatment, levels = c(1, 2))
boots_pairs_w_sup$label <- paste0(substr(boots_pairs_w_sup$focal, 1, 2), "-", substr(boots_pairs_w_sup$competitor, 1, 2))
write.csv(boots_pairs_w_sup, "./output/boots_pairs_w_sup.csv")

# Create coexistence area for plot - min/max fitness difference that permits coexistence
niche_differentiation <- seq(from = -.25, to = 1, by = 0.001)
niche_overlap <- 1-niche_differentiation
fitness_ratio_min <- niche_overlap
fitness_ratio_max <- 1/niche_overlap

df <- data.frame(niche_diff = niche_differentiation,
                 min_fitness_ratio = fitness_ratio_min,
                 max_fitness_ratio = fitness_ratio_max)
head(df)

#Text size 
geom.text.size = 7*(5/14)
geom.text.size2 = 6*(5/14)
label.size = 5*(5/14)
element.size = 5
theme.size = 7

#Coexistence plot: Figure 1
ggplot(boots_pairs_w_sup) + 
  theme_classic()+
  theme(text = element_text( size = theme.size),
        axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 1)),
        axis.text.y = element_text(size = (theme.size - 1)),
        legend.title = element_text(size = theme.size-1),
        legend.text = element_text(size = theme.size-1),
        legend.position = c(.2, .01),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        #legend.margin = margin(.5, .5, .5, .5),
        legend.direction = "horizontal",
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.x = unit(-0.2, "cm"),
        plot.margin=margin(t = 3, r = 0, b = 2, l = 2)) +
  coord_cartesian(xlim = c(-.01, 1), ylim = c(0.25, 10)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  annotate( "text", x = 0.8, y = 0.9, label = "Coexistence: \n \u03c1 < (Kj/Ki) < 1/\u03c1", size = geom.text.size2, fontface = "bold") +
  annotate("text", x = 0.13, y = 6.5, label = "Competitive Exclusion", size = geom.text.size2, fontface = "bold") +
  annotate("text", x = 0.13, y = .55, label = "Competitive Exclusion", size = geom.text.size2, fontface = "bold") +
  geom_point(aes(x = snd, y = fd, color = treat, shape = treat, group = sp_pair), size = 0.75, stroke = 1) +
  scale_color_manual(values=c('1' ="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1'= 15, '2'= 16), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  geom_text_repel(data = boots_pairs_w_sup %>% filter(treat == 2), 
                  aes(x = snd, y = fd, label = label), 
                  size = label.size,
                  box.padding = 0.05,
                  nudge_x = c(0.025, 0.025, 0.01, 0.025, .025, 0,.025, -0.025, -.025, -0.025, -.025,-0.025,0.025,0.025,0.025 ), 
                  nudge_y = c(0.005, 0, 0.002, 0, 0,0.05, 0, -0.02, 0, 0, 0.01, 0,0,-.005,0 )) +
  geom_path(aes(x = snd, y = fd, group = sp_pair), arrow = arrow(length=unit(0.03,"inches"), ends="last", type = "closed"), size = 0.25) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(expand= c(.1,0)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)", size = theme.size) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1,byrow = TRUE), col = guide_legend(nrow = 2) , size = F)+
  NULL
ggsave("./figures/boots_w_sup.pdf", width = 89, height = 55, units = "mm", device = cairo_pdf)

#lambda plot: Figure 2
seed_data$treat <- factor(seed_data$treat, levels = c("W", "D"))

seed_data %>%
  filter(num_comp == 0) %>%
  ggplot(aes(x = focal, y = num_seeds, fill = factor(treat, levels = c("W", "D")))) +
  scale_x_discrete(labels = c( "ACWR" = "AC", "FEMI" = "FE", "HOMU" = "HO", "PLER" = "PL", "SACO" = "SA", "URLI"= "UR")) +
  theme_classic() +
  theme(text = element_text(size = theme.size),
        axis.text.y = element_text(size = (theme.size - 1)),
        legend.position = c(.82, 1),
        legend.justification = c("top"),
        legend.box.just = "top",
        legend.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = theme.size-1),
        legend.text = element_text(size = theme.size-1), 
        legend.key.size = unit(.3, "cm"),
        plot.margin=margin(t = 3, r = 3, b = 0, l = 3)) +
  guides(fill = guide_legend(title = "Treatment:", title.position = "left"), 
         col = guide_legend(nrow = 1), 
         size = "none" )+
  geom_boxplot(width=0.5, lwd=0.25,position=position_dodge(.6), outlier.size=0.25)+
  #geom_point(aes(shape = treat, group = treat), position=position_jitterdodge(dodge.width=0.9), size=2)+
  scale_fill_manual(values=c("W" ="#4E84C4", "D" = "#D16103"), name = "treat", labels = c("Ambient", "Reduced Rain")) +
  geom_text(aes(x = 1, y = 10000, label = "ns"), size =label.size) +
  geom_text(aes(x = 2, y = 5000, label = "*"), size = 8*(5/14)) +
  geom_text(aes(x = 3, y = 3000, label = "ns"), size = label.size) +
  geom_text(aes(x = 4, y = 3000, label = "ns"), size = label.size) +
  geom_text(aes(x = 5, y = 4000, label = "ns"), size = label.size) +
  geom_text(aes(x = 6, y = 1400, label = "*"), size =8*(5/14)) +
  xlab(" ")+
  ylab("Fecundity")+
  labs( color = "treat") +
  scale_y_log10() +
  NULL
ggsave("./figures/real_lambda_log.pdf", width = 89, height = 55, units = "mm")


#TRAITS plot: Figure 4
  

  df_compile<- read_csv("./output/pca_mantel.csv")
#source("./scripts/pca_traits.R")
  #these are from pca_traits.R (run first or source above) - they provide r2 and p values
  vegan :: mantel(dist_fd_d, dist_pc1)
  vegan :: mantel(dist_fd_w, dist_pc1)
  vegan :: mantel(dist_snd_d, dist_pc1)
  vegan :: mantel(dist_snd_w, dist_pc1)
  vegan :: mantel(dist_fd_diffs_abs, dist_pc1)
  vegan :: mantel(dist_snd_diffs, dist_pc1)
  
 fd_plot <- ggplot(df_compile)+
    geom_point(aes(x = pc1, y = fd_w ), color = "#4E84C4", shape = 15, size = 0.75) +
    geom_point(aes(x = pc1, y = fd_d), color = "#D16103", shape = 19, size = 0.75) +
    theme_classic() +
   theme(text = element_text(size = theme.size))+
   annotate("text", x = 0.5, y = 8.6, label = "Reduced Rain", size = geom.text.size, color = "#D16103", hjust = 0) +
   annotate("text", x = 0.5, y = 8, label = bquote(italic(R)^2 ~.("= 0.839")), size = geom.text.size, color = "#D16103", hjust = 0) +
    annotate("text", x = 0.5, y = 7.4, label = "p = 0.003*", size = geom.text.size, color = "#D16103",  hjust = 0) +
   annotate("text", x = 4, y = 2.3, label = "Ambient", size = geom.text.size, color = "#4E84C4",  hjust = 0) +
   annotate("text", x = 4, y = 1.5, label = bquote(italic(R)^2 ~.("= 0.423")), size = geom.text.size, color = "#4E84C4", hjust = 0) +
    annotate("text", x = 4, y = 0.8, label = "p = ns", size = geom.text.size, color = "#4E84C4", hjust = 0) +
    stat_smooth(aes(x = pc1, y = fd_w), color = "#4E84C4", method = "lm", se=F) +
    stat_smooth(aes(x = pc1, y = fd_d), color = "#D16103", method = "lm", se = F, linetype = "dashed") +
    ylab("Fitness differneces") +
    xlab("functional trait distances \n between pairs on pc1 axis")
  
  snd_plot <- ggplot(df_compile)+
    geom_point(aes(x = pc1, y = snd_w ), color = "#4E84C4", shape = 15, size = 0.75) +
    geom_point(aes(x = pc1, y = snd_d), color = "#D16103", shape = 19, size = 0.75) +
    theme_classic() +
    theme(text = element_text(size = theme.size))+
    annotate("text", x = 0.5, y = 1.1, label = "Reduced Rain", size = geom.text.size, color = "#D16103", hjust = 0) +
    annotate("text", x = 0.5, y = 1, label = bquote(italic(R)^2 ~.("= 0.266")), size = geom.text.size, color = "#D16103", hjust = 0) +
    annotate("text", x = 0.5, y = 0.9, label = "p = ns", size = geom.text.size, color = "#D16103",  hjust = 0) +
    annotate("text", x = 4, y = 0.2, label = "Ambient", size = geom.text.size, color = "#4E84C4",  hjust = 0) +
    annotate("text", x = 4, y = 0.1, label = bquote(italic(R)^2 ~.("= 0.031")),  size = geom.text.size, color = "#4E84C4", hjust = 0) +  
    annotate("text", x = 4, y = 0.02, label = "p = ns", size = geom.text.size, color = "#4E84C4",  hjust = 0) +
    stat_smooth(aes(x = pc1, y = snd_w), color = "#4E84C4", method = "lm", se = F) +
    stat_smooth(aes(x = pc1, y = snd_d), color = "#D16103", method = "lm", se = F, linetype = "dashed") +
    ylab("Stabilizing niche differneces") +
    xlab("functional trait distances \n between pairs on pc1 axis") 
  
  snd_diff_plot <- ggplot(df_compile, aes(x = pc1, y = snd_diffs))+
    geom_point( size = 0.75)+
    theme_classic() +
    theme(text = element_text(size = theme.size))+
    geom_text_repel(aes(label = label), size = geom.text.size2,
                    box.padding = 0.1,
                    nudge_y = c(0.01, 0.01, 0, 0.01, -.01, .01 , .01, 0.01, -.01, -.01, .01, .01, -.005, -.001, .01), 
                    nudge_x = c(.1, .1, .05, -.1, .1, .1 , .1, 0.05, -.01, .1, .1, .1, 0.1, 0.1, .1)) +
    annotate("text", x = 0.5, y = .68, label = bquote(italic(R)^2 ~.("= 0.167")), size =geom.text.size, hjust = 0) +
    annotate("text", x = 0.5, y = .62, label = "p = ns", size = geom.text.size, hjust = 0) +
    stat_smooth(method = "lm", color = "black", se = F) +
    ylab("Change in stabilizing niche \n differneces between treatments") +
    xlab("functional trait distances \n between pairs on pc1 axis") 
  
  fd_diff_plot_abs <- ggplot(df_compile, aes(x = pc1, y = fd_diff_abs))+
    geom_point(size = 0.75)+
    theme_classic() +
    theme(text = element_text(size = theme.size))+
    geom_text_repel(aes(label = label), size = geom.text.size2,
                    box.padding = 0.1,
                    nudge_y = c(.1, -.1, -.01, .1, .1, -.1 , .1, 0, .1, -.1, .1, .1, .1, .05, .1), 
                    nudge_x = c(.1, .1, .1, .1, .1, .1 , .1, 0, -.1, .1, .1, -.05, 0.1, 0.1, .1))+
    annotate("text", x = 0.5, y = 5.1, label = bquote(italic(R)^2 ~.("= 0.537")), size = geom.text.size, hjust = 0)+
    annotate("text",x = 0.5, y = 4.55, label = "p = 0.028*", size = geom.text.size, hjust = 0) +
    stat_smooth(method = "lm", color = "black", se= F) +
    ylab("Change in fitness differneces \n between treatments (absolute value)") +
    xlab("functional trait distances \n between pairs on pc1 axis") 
  
diff_plots <- ggarrange(fd_diff_plot_abs, snd_diff_plot, labels = c("a", "b"), font.label = list(size = 8, color = "black", face = "bold"))
fd_snd_plot <- ggarrange(fd_plot, snd_plot, labels = c("c", "d"), font.label = list(size = 8, color = "black", face = "bold"))
ggarrange( diff_plots, fd_snd_plot, nrow = 2)
  
ggsave("./figures/trait_diff.pdf", width = 183, height = 110, units = "mm", device = cairo_pdf)

#Supplemental Coexistence plot: Figure ED1

ggplot(boots_pairs_w_sup) + 
  theme_classic()+
  theme(text = element_text(size = theme.size))+
  coord_cartesian(xlim = c(-.05, 1), ylim = c(0.1, 15)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treat, shape = treat, group = sp_pair), size = 2, stroke = 1) +
  geom_errorbar(aes(x = snd, ymin = (fd_low), ymax = (fd_high))) +
  geom_errorbarh(aes(y = fd, xmin = (snd_low), xmax = (snd_high))) +
  scale_color_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1' = 15, '2' = 16), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)") +
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
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        panel.spacing = unit(0,'lines')) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , scale = "none")+
 facet_wrap(vars(label), nrow = 4, scales ="free")+
  NULL
ggsave("./figures/coexistence_facet.jpg", width = 183, height = 183,units ="mm", device = "jpeg")


#eta_alpha plot: Figure 3

n_alpha_ratios <- read.csv("./output/n_alph_ratio_output.csv")

ggplot(n_alpha_ratios, aes(x = type_change, y =change_value))+
  theme_classic() +
  theme(text = element_text(size = theme.size))+
  geom_boxplot(fill = "light grey") +
  ylab("Difference between treatments") +
  ylim(c(0,1)) +
  xlab (" ") +
  scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) +
  geom_text(aes(x = 1.5, y = 0.9, label = "*"), size = 8) +
  NULL
ggsave("./figures/alpha_eta_ratio_box.pdf", width= 89, height = 70, units = "mm", device = cairo_pdf)


### TABLES------
library(gt)

##Stabilizing niche and fitness differences: Table ED3
pars_boot <- nls_boot_pairs_unique[, c("focal","competitor", "treatment", "snd", "fd", "fd_superior", "coexist", "fd_sup_sp")]
pars_boot$species <- paste0(substr(pars_boot$focal, 1, 2), "-", substr(pars_boot$competitor, 1, 2))

pars_boot$outcome <- ifelse(pars_boot$coexist == 1, "coexist",
                            ifelse(pars_boot$fd_sup_sp == 2, paste0(substr(pars_boot$species, 4, 5), " wins"), 
                                   paste0(substr(pars_boot$species, 1, 2), " wins")))


pair_labels <- pars_boot %>%
  filter(treatment == 1 ) %>%
  filter(fd > 1)
pair_labels <- unique(pair_labels$species)
pars_boot <- pars_boot %>% filter(species %in% pair_labels)
pars_boot <- subset(pars_boot, select =  -c(focal, competitor, coexist, fd_sup_sp, fd_superior))

pars_boot_wide<-pivot_wider(data = pars_boot,
                        names_from = treatment,
                        values_from = c(snd, fd, outcome))
pars_boot_wide <- pars_boot_wide %>% arrange(species)
pars_tab <- gt(pars_boot_wide) 

pars_tab <- pars_boot_wide %>%
  gt(rowname_col ="species") %>%
  tab_stubhead(label = "Species Pair") %>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = c(snd_1, fd_1, outcome_1)) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = c(snd_2, fd_2, outcome_2)) %>%
  cols_move_to_start(
    columns = c(snd_1, fd_1, outcome_1)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
  tab_style(
    style = cell_text(align = "center"), 
    locations = cells_stubhead()) %>%
  cols_width(species~px(100), snd_1~px(80), fd_1~px(80), snd_2~px(80), fd_1~px(80)) %>%
  cols_label(
   snd_1 = html("Stabilizing<br>niche<br>difference"), 
   fd_1 = html("Fitness<br>difference"), 
   outcome_1 = html("Predicted<br>Outcome"), 
   snd_2 = html("Stabilizing<br>niche<br>difference"), 
   fd_2 = html("Fitness<br>difference"), 
   outcome_2 = html("Predicted<br>Outcome")
  ) %>%
  cols_align(align = "center") %>%
  fmt_number(c(snd_1, fd_1, snd_2, fd_2), decimals = 3) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_stub()) %>%
  tab_row_group(
    label = html("Coexist in<br>ambient but not<br>reduced rainfall"), 
    rows = (outcome_1 == "coexist" & outcome_2 != "coexist"), 
    id = "aa") %>%
  tab_row_group(
    label = html("Coexist in<br>reduced rainfall<br>but not ambient"), 
    rows = c(outcome_2 == "coexist" & outcome_1 != "coexist"), 
    id = "bb") %>%
  tab_row_group(
    label = html("Coexist in<br>both"), 
    rows = (outcome_1 == "coexist" & outcome_2 == "coexist"), 
    id = "cc") %>%
  tab_row_group(
    label = html("Coexist in<br>neither"), 
    rows = (outcome_1 != "coexist" & outcome_2 != "coexist"),
    id = "dd") %>%
  row_group_order(
    groups = c("aa", "bb", "cc", "dd")
  ) %>%
  tab_options(row_group.as_column = T)%>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_row_groups())

pars_tab

gtsave(pars_tab, "snd_fd_table.png", "./figures/", expand = T)

##Lambda Table : Table ED1

lambda <- seed_data %>%
  filter(num_comp == 0) 

lambda_tab_df<- lambda %>%
  dplyr::group_by(focal, treat) %>%
  dplyr::mutate(count=n()) %>%
  dplyr::summarize(mean  = mean(num_seeds, na.rm = TRUE), se = sd(num_seeds, na.rm = TRUE)/sqrt(length(num_seeds)), count = mean(count))

lambda_tab_wide <- pivot_wider(data = lambda_tab_df, 
                               names_from = treat, 
                               values_from = c(mean, se, count))
# Add p-values from glmer fecundity ~ species*treatment + plot -- found in lambda_plants.R
lambda_tab_wide$p_value <- c(0.3597, "0.0006*", 0.6869, 0.0670, 0.4452,"<0.0001*" )
#Edit other columns
lambda_tab_wide$species_code<- c("AC", "FE", "HO", "PL", "SA", "UR")
lambda_tab_wide$sci_name <- c("Acmispon wrangelianus", "Festuca microstachys", "Hordeum murinum", "Plantago erecta", "Salvia columbariae", "Uropappus lindleyi")
lambda_tab_wide$family <- c("Fabaceae", "Poaceae", "Poaceae", "Plantaginaceae", "Lamiaceae", "Asteraceae")
#need to get rid of the group_by legacy
lambda_tab_wide <- as.data.frame(lambda_tab_wide)
lambda_tab_wide$focal <- NULL

#Table
lambda_tab <- lambda_tab_wide %>%
  gt(rowname_col = "family") %>%
  tab_stubhead(label = "Family") %>%
  tab_style( 
    style = cell_text(align = "left"), 
    locations = cells_stub())%>%
  tab_style( 
    style = cell_text(align = "left", style = "italic"), 
    locations = cells_body(columns = "sci_name") )%>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = c(mean_W, se_W, count_W)) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = c(mean_D, se_D, count_D)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
   cols_move_to_start(
    columns = c(sci_name, species_code)) %>%
  cols_move(
    columns = c(mean_W, se_W, count_W), 
    after = species_code) %>%
  cols_label(
    sci_name = "Species",
    species_code = "Code",
    mean_W = "Mean Fecundity", 
    se_W = "Standard error", 
    count_W = "n",
    mean_D = "Mean Fecundity", 
    se_D = "Standard error", 
    count_D = "n",
    p_value = "p-value") %>%
  cols_align(align = "center") %>%
  cols_width(family~px(120), sci_name~px(100), species_code~px(70), count_D~px(50), count_W~px(50), everything()~px(90)) %>%
  fmt_number(c(mean_W, se_W, mean_D, se_D), decimals = 1)
  
lambda_tab

gtsave(lambda_tab, "lambda_tab.png", "./figures/")

#Gravimetric Water Content Table: Table ED2

gwc <- read.csv("./data/soil_gwc.csv")
gwc$GWC <- ((gwc$wet_weight - gwc$dry_weight)/gwc$dry_weight)

#Get the df ready for table 
gwc<- gwc %>% 
  dplyr::group_by(doy, treatment) %>%
  dplyr::summarize(mean = mean(GWC))

gwc_wide <- pivot_wider(data = gwc, 
                               names_from = treatment, 
                               values_from = mean)
gwc_wide$date <- c("March 27, 2019", "April 21, 2019", "May 17, 2019")
gwc_wide$p_value <- c(0.0005, 0.00004, 0.017) #from soil_gwc.R
gwc_wide <- as.data.frame(gwc_wide)
gwc_wide$doy <- NULL

#Make Table
gwc_tab <- gwc_wide %>%
  gt(rowname_col ="date") %>%
  tab_stubhead(label = "Date") %>%
  cols_label(
    W = "Mean GWC", 
    D = "Mean GWC", 
    p_value = "p-value") %>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = W) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = D) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_stub()) %>%
  cols_align(align = "center") %>%
  cols_move_to_start( W) %>%
  fmt_number(c(W, D), decimals = 3) %>%
  fmt_number(p_value, decimals = 5)
gwc_tab
gtsave(gwc_tab, "gwc_table.png", "./figures/")

##Trait table: table ED4

trait_df <- read.csv("./data/trait_list.csv")

trait_tab <- trait_df %>%
  gt() %>%
  tab_style(cell_text(align = "left"), 
            locations = cells_stub()) %>%
  tab_style(cell_text(align = "center"), 
            locations = cells_body(columns = "Units")) %>%
  tab_row_group(
    label = "Leaf", 
    rows = (Organ == "leaf"), 
    id = "aa") %>%
  tab_row_group(
    label = "Seed", 
    rows = (Organ == "seed"), 
    id = "bb") %>%
  tab_row_group(
    label = "Root", 
    rows = (Organ == "root"), 
    id = "cc") %>%
  tab_row_group(
    label = html("Whole<br>plant"), 
    rows = (Organ == "whole_plant"), 
    id = "dd") %>%
  row_group_order(
    groups = c("aa", "cc", "dd", "bb")) %>%
  tab_options(row_group.as_column = T)%>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_row_groups()) %>%
  cols_hide(columns = "Organ") %>%
  tab_style(
    style = cell_borders(
      sides = c("left"),
      color = "grey",
      weight = px(1),
      style = "solid"
    ),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )) %>%
  text_transform(
    locations = cells_body(
      columns = "Units",
      rows = 1
    ),
    fn = function(x){
      paste0("cm", "&sup2")
    }) %>%
  text_transform(
    locations = cells_body(
      columns = "Units",
      rows = 2
    ),
    fn = function(x){
      paste0( "<sup>cm&sup2</sup>&frasl;<sub>g</sub>")
    }) %>%
  text_transform(
    locations = cells_body(
      columns = "Units",
      rows = c(3, 4)
    ),
    fn = function(x){
      paste0( "<sup>mg</sup>&frasl;<sub>g</sub>")
    }) %>%
  text_transform(
    locations = cells_body(
      columns = "Units",
      rows = 7
    ),
    fn = function(x){
      paste0( "<sup>m</sup>&frasl;<sub>g</sub>")
    }) %>%
  text_transform(
    locations = cells_body(
      columns = "Units",
      rows = 11
    ),
    fn = function(x){
      paste0("d<sub>13</sub>C")
    })

trait_tab

gtsave(trait_tab, "trait_description_table.png", "./figures/")

