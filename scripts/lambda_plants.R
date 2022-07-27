library(lme4)
library(ggfortify)
library(emmeans)
library(tidyverse)

seed_data <- read.csv("./data/drought_seed_production_data.csv")
seed_data$treat <- factor(seed_data$treat, levels = c("W", "D"))
lamda_plants <- seed_data %>% filter(num_comp < 1)
lamda_plants$num_seeds <- floor(lamda_plants$num_seeds)
out1 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "ACWR"))
out2 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "FEMI"))
out3 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "HOMU"))
out4 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "PLER"))
out5 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "SACO"))
out6 <- glmer(num_seeds ~ treat + (1|plot), family = poisson(link = "log"), data = lamda_plants %>% filter(focal == "URLI"))
summary(out1);summary(out2);summary(out3);summary(out4);summary(out5);summary(out6)
plot(out1)
qqnorm(resid(out3))
qqline(resid(out3))



#What's the right distribution to use?

lamda_plants %>% filter(focal == "PLER") %>%
  ggplot() +
  geom_histogram(aes(x = num_seeds), binwidth = 25)
with(lamda_plants, mean(num_seeds[focal== "PLER"]))

seed_data %>%
  filter(num_comp == 0) %>%
  ggplot(aes(x = focal, y = num_seeds, fill = factor(treat))) +
  geom_boxplot() +
  scale_fill_manual(values=c("W" ="#00BFC4", "D" = "#F8766D"), name = "treat", labels = c("W", "D")) +
  #ylim(0, 7000) +
  ggtitle("Plants grown without competitors")+
  xlab(" ")+
  ylab("fecundity")+
  theme(axis.title = element_text(size = 18), title = element_text(size = 20))+
  #stat_summary(fun=mean, geom="point", shape=18, size=4, position=position_dodge(.75))+
  NULL
ggsave("./figures/real_lambda_log.pdf", width = 11, height = 6)

seed_data %>%
  filter(num_comp == 0) %>%
  ggplot(aes(x = focal, y = num_seeds, fill = factor(treat))) +
  geom_boxplot() +
  scale_fill_manual(values=c("W" ="#00BFC4", "D" = "#F8766D"), name = "treat", labels = c("W", "D")) +
  ylim(0, 7000) +
  ggtitle("Plants grown without competitors")+
  xlab(" ")+
  ylab("fecundity")+
  theme(axis.title = element_text(size = 18), title = element_text(size = 20)) +
  #stat_summary(fun=mean, geom="point", shape=18, size=4, position=position_dodge(.75))+
  NULL
ggsave("./figures/real_lambda_no_out.pdf", width = 11, height = 6)


#Stats Consultant help: https://stats.idre.ucla.edu/r/seminars/interactions-r/#s5
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

seed_data$treat <- factor(seed_data$treat, levels = c("W", "D"))

seed_data %>%
  filter(num_comp == 0) %>%
  ggplot(aes(x = focal, y = log(num_seeds), fill = factor(treat))) +
  theme_classic(base_size = 20) +
  theme( plot.title = element_text(size = 20), 
         plot.subtitle = element_text(size = 15), 
         axis.title = element_text(size = 18)) +
  geom_boxplot() +
  stat_summary()  +
  scale_fill_manual(values=c("W" ="lightblue2", "D" = "#D55E00"), name = "treat", labels = c("W", "D")) +
  geom_text(aes(x = 1, y = 2.25, label = "p=0.3597"), size = 5) +
  geom_text(aes(x = 2, y = 3.75, label = "p=0.0006*"), size = 5) +
  geom_text(aes(x = 3, y = 3.75, label = "p=0.6869"), size = 5) +
  geom_text(aes(x = 4, y = 3.75, label = "p=0.0670"), size = 5) +
  geom_text(aes(x = 5, y = 4, label = "p = 0.4452"), size = 5) +
  geom_text(aes(x = 6, y = 3.25, label = "p = < 0.0001*"), size = 5) +
  labs(title = "Plants grown without competitors", 
       subtitle = "GLMER Model: fecundity ~ species*treatment + plot") +
  xlab(" ")+
  ylab("ln(fecundity)")+
  NULL

find_se <- function(trmt, species){
aa <- mean(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] )
ss <- sd(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] ) 
nn <- length(seed_data$num_seeds[seed_data$num_comp == 0 & seed_data$focal == species & seed_data$treat == trmt ] ) 
se <- ss/sqrt(nn)
print (ss)
print(aa)
print(se)}

find_se("W", "FEMI")
find_se("W", "URLI")
