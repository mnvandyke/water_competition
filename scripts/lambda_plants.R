### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### this script tests differences in seed production between treatments in plots with zero competitors. 

library(lme4)
library(ggfortify)
library(emmeans)
library(tidyverse)

seed_data <- read.csv("./data/drought_seed_production_data.csv")
seed_data$treat <- factor(seed_data$treat, levels = c("W", "D"))
lamda_plants <- seed_data %>% filter(num_comp < 1)
lamda_plants$num_seeds <- floor(lamda_plants$num_seeds)

#Using UCLA Stats Consultant Strategy: https://stats.idre.ucla.edu/r/seminars/interactions-r/#s5
lamda_plants$focal <- as_factor(lamda_plants$focal)
lamda_plants$focal <- relevel(lamda_plants$focal, ref = "ACWR")
lamda_plants$treat <- as_factor(lamda_plants$treat)
lamda_plants$treat <- relevel(lamda_plants$treat, ref = "W")

focal_trt <- glmer(num_seeds~focal*treat + (1|plot), family = poisson(link = "log"), data = lamda_plants)
summary(focal_trt) 

em_focal_trt <- emmeans(focal_trt, ~focal*treat )
summary(em_focal_trt)
contrast(em_focal_trt, "pairwise", by = "focal", adjust = "tukey")
emmip(focal_trt, focal ~ treat, CIs = TRUE)


find_se <- function(trmt, species){
aa <- mean(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] )
ss <- sd(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] ) 
nn <- length(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] ) 
se <- ss/sqrt(nn)
print (ss)
print(aa)
print(se)}

find_se("W", "FEMI")
find_se("D", "FEMI")

#Check that the estimated lambdas using all plots is similar

final_output <- read.csv("./output/final_output_nls_boot_1000.csv")
final_output$treatment <- factor(final_output$treatment, levels = c(1, 2))

final_output %>%
  ggplot(aes(x = focal, y = lambda, fill = factor(treatment))) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values=c("1" ="#4E84C4", "2" = "#D16103"), name = "treatment", labels = c("Ambient", "Reduced Rain")) 
  


