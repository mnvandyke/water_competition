### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### this script makes the structural analysis tables (Extended data tables 5,6,7)

library(gt)
library(tidyverse)
library(dplyr)


struc <- read.csv("./output/structural_analysis_output.csv") %>% arrange(species)
struc$X <- NULL
struc$coexist <- ifelse(struc$outcome == 1, "yes", "no")
struc$outcome <- NULL

##Triplets Table : table ED6

trip <- struc %>%
  filter(no_sp == 3)
trip$match <- NULL
trip$no_sp <- NULL

trip_wide<-pivot_wider(data = trip,
            names_from = treatment,
            values_from = c(theta, omega, coexist)
            )


trip_tab <- gt(trip_wide) 
trip_tab <- trip_wide %>%
  gt(rowname_col ="species") %>%
  tab_stubhead(label = "Species") %>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = c(omega_1, theta_1, coexist_1)) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = c(omega_2, theta_2, coexist_2)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
  cols_label(
    omega_1 = html("&#937"), 
    theta_1 = html("&#952"), 
    coexist_1 = html("Predicted <br> Coexistence?"), 
    omega_2 = html("&#937"), 
    theta_2 = html("&#952"), 
    coexist_2 = html("Predicted <br> Coexistence?")
  ) %>%
  cols_align(align = "center") %>%
  #cols_width(species~px(110)) %>%
  fmt_number(c(omega_1, omega_2), decimals = 3) %>%
  fmt_number(c(theta_1, theta_2), decimals = 2) %>%
  tab_style(
    style = cell_text(align = "center", indent = px(100)), 
    locations = cells_stubhead()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_stub())%>%
  tab_row_group(
    label = html("Coexist in <br> ambient but not <br> reduced rainfall"), 
    rows = (coexist_1 == "yes" & coexist_2 == "no"), 
    id = "aa") %>%
  tab_row_group(
    label = html("Coexist in <br> reduced rainfall <br> but not ambient"), 
    rows = (coexist_1 == "no" & coexist_2 == "yes"), 
    id = "bb") %>%
  tab_row_group(
    label = html("Coexist in <br> neither"), 
    rows = (coexist_1 == "no" & coexist_2 == "no"), 
    id = "cc") %>%
  row_group_order(
    groups = c("aa", "bb", "cc")) %>%
  tab_options(row_group.as_column = T)%>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_row_groups())
trip_tab
gtsave(trip_tab, "structural_table_triplets.png", "./figures/")


## Pairs Table: Table ED5
pairs <- struc %>%
  filter(no_sp == 2)
pairs$no_sp <- NULL
pairs <- pairs %>% filter(species %in% pair_labels) ## pairs_labels is from final_figures.R

pairs_wide<-pivot_wider(data = pairs,
                       names_from = treatment,
                       values_from = c(theta, omega, coexist)
)

pairs_tab <- gt(pairs_wide) 
pairs_tab <- pairs_wide %>%
  gt(rowname_col ="species") %>%
  tab_stubhead(label = "Species") %>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = c(omega_1, theta_1, coexist_1)) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = c(omega_2, theta_2, coexist_2)) %>%
  cols_move_to_start(
    columns = c(omega_1, theta_1, coexist_1)) %>%
  #cols_width(species~px(110), omega_1~px(100), omega_2~px(100), theta_1~px(100), theta_2~px(100)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
  cols_label(
    omega_1 = html("&#937"), 
    theta_1 = html("&#952"), 
    coexist_1 = html("Predicted <br> Coexistence?"), 
    omega_2 = html("&#937"), 
    theta_2 = html("&#952"), 
    coexist_2 = html("Predicted <br> Coexistence?")) %>%
  cols_align(align = "center") %>%
  fmt_number(c(omega_1, omega_2), decimals = 3) %>%
  fmt_number(c(theta_1, theta_2), decimals = 2) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_stub()) %>%
  tab_style(
    style = cell_text(align = "center", indent = px(125)), 
    locations = cells_stubhead()) %>%
  tab_row_group(
    label = html("Coexist in <br> ambient but not <br> reduced rainfall"), 
    rows = (coexist_1 == "yes" & coexist_2 == "no"), 
    id = "aa") %>%
  tab_row_group(
    label =html("Coexist in<br>reduced rainfall<br>but not ambient"), 
    rows = (coexist_1 == "no" & coexist_2 == "yes"), 
    id = "bb") %>%
  tab_row_group(
    label = html("Coexist in <br> both"), 
    rows = (coexist_1 == "yes" & coexist_2 == "yes"), 
    id = "cc") %>%
  tab_row_group(
    label = html("Coexist in <br> neither"), 
    rows = (coexist_1 == "no" & coexist_2 == "no"),
    id = "dd") %>%
  row_group_order(
    groups = c("aa", "bb", "cc", "dd")) %>%
  tab_options(row_group.as_column = T) %>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_row_groups()
  )

pairs_tab

gtsave(pairs_tab, "structural_table_pairs.png", "./figures/")

##Quads, Quints, Sexts Table: Table ED7

quads <- struc %>%
  filter(no_sp >= 4)
quads_wide<-pivot_wider(data = quads,
                        names_from = treatment,
                        values_from = c(theta, omega, coexist)
)

#quads_tab <- gt(quads_wide) 
quads_tab <- quads_wide %>%
  gt(rowname_col ="species") %>%
  tab_stubhead(label = "Species") %>%
  tab_spanner(
    label = "Ambient Rainfall", 
    columns = c(omega_1, theta_1, coexist_1)) %>%
  tab_spanner(
    label = "Reduced Rain", 
    columns = c(omega_2, theta_2, coexist_2)) %>%
  cols_move_to_start(
    columns = c(omega_1, theta_1, coexist_1)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_spanners(), cells_stubhead())) %>%
  cols_label(
    omega_1 = html("&#937"), 
    theta_1 = html("&#952"), 
    coexist_1 = html("Predicted <br> Coexistence?"), 
    omega_2 = html("&#937"), 
    theta_2 = html("&#952"), 
    coexist_2 = html("Predicted <br> Coexistence?")
  ) %>%
  cols_align(align = "center") %>%
  #cols_width( omega_1~px(90), omega_2~px(90), 
             #theta_1~px(90), theta_2~px(90), coexist_1~px(90), coexist_2~px(60)) %>%
  fmt_number(c(omega_1, omega_2), decimals = 4) %>%
  fmt_number(c(theta_1, theta_2), decimals = 2) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_stub()) %>%
  tab_style(
    style = cell_text(align = "center"), 
    locations = cells_stubhead()) %>%
  tab_style(
    style = cell_text(align  = 'left'), 
    locations = cells_row_groups()
  ) %>%
  tab_row_group(
    label = "Quadruplets", 
    rows = no_sp == 4, 
    id = "quads") %>%
  tab_row_group(
    label = "Quintuplets", 
    rows = no_sp == 5, 
    id = "quints") %>%
  tab_row_group(
    label = "Sextuplet", 
    rows = no_sp == 6, 
    id = "sext") %>%
  row_group_order(c("quads", "quints", "sext")) %>%
  cols_hide(no_sp) %>%
  tab_options(row_group.as_column = T)%>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_row_groups())
quads_tab  

gtsave(quads_tab, "structural_table_quads.png", "./figures/")
