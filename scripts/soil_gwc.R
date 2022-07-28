### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### this script tests for differences in gravimetric water content between treatments 
### at three time points during the growing season

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
options(scipen = 5)

soil <- read.csv("./data/soil_gwc.csv")


soil$gwc <-((soil$wet_weight - soil$dry_weight)/soil$dry_weight)

ggplot(soil, aes(x = treatment, y = gwc))+
  geom_boxplot()

april <- soil %>%
  filter(date == "4/21/19")

may <- soil %>%
  filter(date == "5/17/19")

march <- soil %>%
  filter(date == "3/27/19")

t.test(soil$gwc ~ soil$treatment)
t.test(april$gwc ~ april$treatment)
t.test(may$gwc ~ may$treatment)
t.test(march$gwc ~ march$treatment)

soil %>% 
  group_by(treatment, date) %>%
  get_summary_stats(gwc, type = "mean_sd")

ggplot(soil, aes(x= date, y = gwc, color = treatment))+
        geom_boxplot()

soil%>%
  group_by(treatment, date) %>%
  identify_outliers(gwc)

ggqqplot(soil, "gwc", ggtheme = theme_bw()) +
  facet_grid(date ~treatment)

test<-aov(gwc ~treatment, data = soil)
summary(test)

