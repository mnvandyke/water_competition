### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 26 July 2022

### This script gets average germination rates for each species 
### by dividing the number of counted germinates by the number of seeds planted

germ <-read.csv("./data/germination.csv")

germ$germ_perc <- 0

germ$germ_perc <- ifelse(germ$focal == "HOMU", germ$germinates/27, ## These numbers are what was planted in an attempt to
                           ifelse(germ$focal == "FEMI", germ$germinates/33,  ## get 25 viable seeds. They are based on viability
                                  ifelse(germ$focal == "PLER", germ$germinates/39, ## tests and are from drought_field_setup.xlsx
                                         ifelse(germ$focal == "SACO", germ$germinates/42, 
                                                ifelse(germ$focal == "URLI", germ$germinates/35, 
                                                       ifelse(germ$focal == "ACWR", germ$germinates/30, 0))))))
                                                                                                                                                 



av_germ <- germ %>%
  group_by(focal) %>%
  summarise(mean_germ = mean(germ_perc, na.rm = T))

av_germ
