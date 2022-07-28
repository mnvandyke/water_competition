### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### This script provides code for the invasion growth rate analysis and output for figure 3


library(dplyr)
library(tidyr)
library(ggplot2)

seed_data <- read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr <- ifelse(seed_data$treat == "W", 1, 2)
seed_data$treatment <- seed_data$Tr
#run nls model script first or upload output
nls_boot_pairs <- read.csv("./output/nls_boot_pairs_1000_full_model.csv")
final_output <- read.csv("./output/final_output_nls_boot_1000.csv")

mean(with(final_output, ni[focal == "HOMU" & treatment == 1]), na.rm = T)


igr_ratio <- function(foc, comp, trt) { #invasion growth rate ratios
  
  nj <- mean(with(final_output, ni[focal == comp & treatment == trt]), na.rm = T)
  ni <- mean(with(final_output, ni[focal == foc & treatment == trt]), na.rm = T)
  ajj <- with(nls_boot_pairs, alpha[focal == comp & competitor == comp & treatment == trt])
  aij <- with(nls_boot_pairs, alpha[focal == foc & competitor == comp & treatment == trt])
  n_ratio = (ni-1)/(nj-1)
  a_ratio = ajj/aij
  return(c( n_ratio, a_ratio, n_ratio*a_ratio))
}
igr_ratio("URLI", "HOMU", 1)

comp_labels <- sort( na.omit( unique(seed_data$background) ))

spp_list <- sort(na.omit( unique(seed_data$focal)))
treat_list <- sort( na.omit( unique(seed_data$treatment)))
spp_treat_comp_combos <- expand.grid(species = spp_list, competitor = comp_labels, treatment = treat_list)

spp_treat_comp_combos$n_ratio <- 0
spp_treat_comp_combos$a_ratio <- 0
spp_treat_comp_combos$product <- 0

for(i in 1:nrow(spp_treat_comp_combos)) {
  ii <- spp_treat_comp_combos[i, "species"] %>% unlist
  jj <- spp_treat_comp_combos[i, "competitor"] %>% unlist
  tt <- spp_treat_comp_combos[i, "treatment"] %>% unlist
  
 spp_treat_comp_combos[i, "n_ratio"] <- igr_ratio(ii, jj, tt)[1]
 spp_treat_comp_combos[i, "a_ratio"] <- igr_ratio(ii, jj, tt)[2]
 spp_treat_comp_combos[i, "product"] <- igr_ratio(ii, jj, tt)[3]
 
}
head(spp_treat_comp_combos)
spp_treat_comp_combos$larger <- ifelse(spp_treat_comp_combos$a_ratio > spp_treat_comp_combos$n_ratio, "a", "n")
spp_treat_comp_combos$invade <- ifelse(spp_treat_comp_combos$product>1, "yes", "no")
head(spp_treat_comp_combos)

invaders <- spp_treat_comp_combos %>%
  filter(invade == "yes")

#Which ratio changes more in invasion growth rate inequality?
igr_change <- function(foc, comp) {
  
  nj_d <- mean(with(final_output, ni[focal == comp & treatment == 2]), na.rm = T)
  ni_d <- mean(with(final_output, ni[focal == foc & treatment == 2]), na.rm = T)
  ajj_d <- with(nls_boot_pairs, alpha[focal == comp & competitor == comp & treatment == 2])
  aij_d <- with(nls_boot_pairs, alpha[focal == foc & competitor == comp & treatment == 2])
  n_ratio_d = log10((ni_d-1)/(nj_d-1))
  a_ratio_d = log10(ajj_d/aij_d)
 
  nj_w <- mean(with(final_output, ni[focal == comp & treatment == 1]), na.rm = T)
  ni_w <- mean(with(final_output, ni[focal == foc & treatment == 1]), na.rm = T)
  ajj_w <- with(nls_boot_pairs, alpha[focal == comp & competitor == comp & treatment == 1])
  aij_w <- with(nls_boot_pairs, alpha[focal == foc & competitor == comp & treatment == 1])
  n_ratio_w = log10((ni_w-1)/(nj_w-1))
  a_ratio_w = log10(ajj_w/aij_w)
  
  nc<-abs(n_ratio_w - n_ratio_d)
  ac<-abs(a_ratio_w - a_ratio_d)
   return(c(nc, ac))
}
igr_change("HOMU", "URLI")

#Create data frame
spp_list <- sort(na.omit( unique(seed_data$focal)))
comp_labels <- sort( na.omit( unique(seed_data$background) ))
spp_comp_combos <- expand.grid(species = spp_list, competitor = comp_labels)

spp_comp_combos$n_change <- 0
spp_comp_combos$a_change <- 0


for(i in 1:nrow(spp_comp_combos)) {
  ii <- spp_comp_combos[i, "species"] %>% unlist
  jj <- spp_comp_combos[i, "competitor"] %>% unlist
 
  
  spp_comp_combos[i, "n_change"] <- igr_change(ii, jj)[1]
  spp_comp_combos[i, "a_change"] <- igr_change(ii, jj)[2]

  
}

spp_comp_combos$larger <- ifelse(abs(spp_comp_combos$a_change)> abs(spp_comp_combos$n_change), "a", 'n')

spp_comp_combos <- spp_comp_combos %>%
  filter(a_change != n_change)
total_pairs <- c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", "URLI_ACWR", 
"HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
"SACO_HOMU", "URLI_HOMU","SACO_PLER", "URLI_PLER", "URLI_SACO")
spp_comp_combos$sp_pair <- paste(spp_comp_combos$species, spp_comp_combos$competitor, sep = "_")

spp_comp_combos <-spp_comp_combos %>%
  filter(sp_pair %in% total_pairs)

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)

pairs_coexist_change <-c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", 
                  "URLI_ACWR",  "HOMU_FEMI", "PLER_FEMI", "URLI_FEMI",
                  "SACO_PLER", "URLI_SACO")

spp_comp_combos_coexist <-spp_comp_combos %>%
  filter(sp_pair %in% pairs_coexist_change)

t.test(spp_comp_combos_coexist$n_change, spp_comp_combos_coexist$a_change, paired = T)


spp_comp_combos_coexist_pivot <- pivot_longer(spp_comp_combos_coexist, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")
  ggplot(spp_comp_combos_coexist_pivot, aes(x = type_change, y =change_value))+
    theme_classic(base_size = 20) +
    theme(text = element_text(size = 18))+
    geom_boxplot() +
    ylab("Difference between treatments (log10)") +
    xlab (" ") +
    scale_x_discrete(labels=c("a_change"="Competitive coefficients \n \u03B1 ratio", "n_change" = "Innate fecundity \n \u03B7 ratio")) +
    geom_text(aes(x = 1.5, y = 1, label = "*"), size = geom.text.size) +
    annotate("text", x = 2, y = 1.5, label = paste("* = p < 0.05"), size = geom.text.size) +
  NULL
  
  t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)
  
  spp_comp_combos_pivot <- pivot_longer(spp_comp_combos, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")
  
  
  ggplot(spp_comp_combos_pivot, aes(x = type_change, y = change_value)) +
    theme_classic(base_size = 20) +
    geom_point(aes(group = sp_pair))+
    geom_path(aes(group = sp_pair)) +
   ylab("Difference between treatments") +
    xlab (" ") +
    scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential"))
  #ggsave("./figures/alpha_eta_ratio_point.pdf")
  
  
  ggplot(spp_comp_combos_pivot, aes(x = type_change, y =change_value))+
    theme_classic(base_size = 20) +
    theme(text = element_text(size = 18))+
    geom_boxplot(fill = "light grey") +
    ylab("Difference between treatments") +
    xlab (" ") +
    #scale_y_log10() +
    scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) +
    geom_text(aes(x = 1.5, y = 1, label = "*"), size = geom.text.size) +
    NULL
ggsave("./figures/alpha_eta_ratio_box.pdf")
  

write.csv(spp_comp_combos_pivot,  "./output/n_alph_ratio_output.csv")
