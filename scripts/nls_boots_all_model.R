### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

###This script fits the annual plant model on 1000 bootstraps of the data 
###and then calculates stabilizing niche and fitness differences for each of those 1000 fits

library(tidyverse)
library(nlstools)
library(dplyr)
library(rsample)
library(broomExtra)
library(purrr)
library(ggplot2)
library(ggrepel)

options(scipen = 5)
set.seed(3)

seed_data <- read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr <- ifelse(seed_data$treat == "W", 1, 2)
seed_data <- seed_data %>% replace_na(list(background = "ACWR")) # putting acwr instead of NA for the lambdas
seed_data$N_acwr <- ifelse(seed_data$background == "ACWR", seed_data$num_comp, 0)
seed_data$N_femi <- ifelse(seed_data$background == "FEMI", seed_data$num_comp, 0)
seed_data$N_homu <- ifelse(seed_data$background == "HOMU", seed_data$num_comp, 0)
seed_data$N_pler <- ifelse(seed_data$background == "PLER", seed_data$num_comp, 0)
seed_data$N_saco <- ifelse(seed_data$background == "SACO", seed_data$num_comp, 0)
seed_data$N_urli <- ifelse(seed_data$background == "URLI", seed_data$num_comp, 0)

homu <- seed_data %>%
  filter(focal == "HOMU")


fit_test<- nls(log(num_seeds)~log(lambda[Tr]/(1+a_ACWR[Tr]*N_acwr+a_FEMI[Tr]*N_femi+a_HOMU[Tr]*N_homu+a_PLER[Tr]*N_pler+a_SACO[Tr]*N_saco+a_URLI[Tr]*N_urli)),
    data=homu, start=list('lambda'= c(100,100), a_ACWR=c(0.1, 0.1), a_FEMI=c(0.1, 0.1), a_HOMU=c(0.1, 0.1), a_PLER=c(0.1, 0.1), a_SACO=c(0.1, 0.1),a_URLI=c(0.1, 0.1)),
    lower = c( 1, 1, rep(.001, 12) ),# setting alpha to be greater than .001  
    upper = c(10000, 10000, rep(2, 12)), 
    control = list(maxiter = 100000),
    algorithm = 'port') 

summary(fit_test)

# function to find lambda and alpha parameter estimates for bootstrapped data, the split is a component of the bootstrap 

fit_nls_on_bootstrap <- function(split) {
  tryCatch({
    nls(log(num_seeds)~log(lambda[Tr]/(1+a_ACWR[Tr]*N_acwr+a_FEMI[Tr]*N_femi+a_HOMU[Tr]*N_homu+a_PLER[Tr]*N_pler+a_SACO[Tr]*N_saco+a_URLI[Tr]*N_urli)),
        analysis(split), 
        start=list('lambda'= c(100,100), a_ACWR=c(0.1, 0.1), a_FEMI=c(0.1, 0.1), a_HOMU=c(0.1, 0.1), a_PLER=c(0.1, 0.1), a_SACO=c(0.1, 0.1),a_URLI=c(0.1, 0.1)),
        lower = c( 1, 1, rep(.001, 12) ),# setting alpha to be greater than .001  
        upper = c(10000, 10000, rep(2, 12)), 
        control = list(maxiter = 1000000),
        algorithm = 'port')},  
    error = function(e) { #Make it keep trying if it fails to converge at first
      print("failed to converge")
    })
}

#Create a data frame with all combinations of treatment and species 
spp_list <- sort(na.omit( unique(seed_data$focal)))
treat_list <- sort( na.omit( unique(seed_data$Tr)))
spp_combos <- expand.grid(species = spp_list)

# loop through the species list, fit the model, return the dataframe 
out <- list() #list for for loop output

for( i in 1:nrow( spp_combos )){ 
  temp_data <- seed_data %>%
    filter(focal == spp_combos[i,1])
  
  #print((temp_data)) 
  
  boots <- bootstraps(temp_data, times = 1300) ## Some boot straps will fail to converge so we do more than necessary
  
  
  boot_models <- boots %>%
    mutate(model = map(splits, fit_nls_on_bootstrap),
           coef_info = map(model, tidy))
  
  boot_coefs <- boot_models %>%
    unnest(coef_info)
  
  boot_coefs$focal <- spp_combos[i, 1]
  
  #print(head(boot_coefs))
  
  out[[i]] <- boot_coefs
  
}


test <- do.call(rbind.data.frame, out)

test<- test %>% 
  select( id, term , estimate, focal) %>% 
  spread(term, estimate) %>% 
  group_by( id)

test<-test%>% pivot_longer(
  cols = starts_with("a"), 
  names_to = c("competitor","treatment"), 
  names_prefix = "a_", 
  names_sep = 4,
  values_to ="alpha"
)

test$lambda<-ifelse(test$treatment == 1, test$lambda1, test$lambda2)
test$lambda1<-NULL
test$lambda2<-NULL

#Now need to throw out the bootstraps that didnt converge and take the first 1000 that did
bootstrap_ids <- unique(test$id)

sp_names <- c("ACWR","FEMI", "HOMU",  "PLER", "SACO",  "URLI")
unique(test$focal[test$id == "Bootstrap0010"]) == sp_names 

vec<-vector()
for( i in bootstrap_ids){
  tt<-ifelse(setequal(sp_names, unique(test$focal[test$id == i])), 0, i)
  vec[i]<- tt
  
}

#vec
false_ids <- vec[vec!=0]
false_ids
test <- subset(test, !(test$id %in% false_ids))
test <-head(test, 72000)
nls_boots <- test


#Check: How many alphas are 0.001?
nrow(nls_boots[nls_boots$alpha == 0.001,])/72000

#Make a data frame for the stabilizing niche and fitness differences calculated from each bootstrap's alphas
boot_list <- sort( na.omit( unique(nls_boots$id)))
spp_treat_boot_combos <- expand.grid(focal = spp_list, treatment = treat_list, competitor = spp_list, id = boot_list)

#Stabilizing Niche Difference -----

stabilizing_niche_diff <- function(df, species1, species2, treat, boot_id) {
  
  
  aij <- with(df, alpha[focal == species1 & competitor == species2 & treatment== treat & id == boot_id])
  #print(paste("aij: ",aij))
  
  aji <- with(df, alpha[focal == species2 & competitor == species1 & treatment== treat & id == boot_id])
  #print(paste("aji: ",aji))
  
  ajj <- with(df, alpha[focal == species2 & competitor == species2 & treatment== treat &  id == boot_id])
  #print(paste("ajj: ",ajj))
  
  aii <- with(df, alpha[focal == species1 & competitor == species1 & treatment== treat & id == boot_id])
  #print(paste("aii: ",aii))
  
  snd <- (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

#test
stabilizing_niche_diff(nls_boots, "ACWR", "URLI", 1, "Bootstrap0010")

#Add snd column to data frame
spp_treat_boot_combos$snd <- 0


for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_boot_combos[i, "competitor"] %>% unlist
  trt <- spp_treat_boot_combos[i , "treatment"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  snd <- stabilizing_niche_diff(nls_boots, sp1, sp2, trt, boot_id)
  spp_treat_boot_combos[i, "snd"] <- snd
  
}

#Getting eta (ηi) term -------
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 

s_g_data <- read.csv("./data/s_g_data.csv") #seed survival and germination data
nls_boots <- merge(nls_boots, s_g_data, by = "focal")
nls_boots$X <- NULL

#ηi equation function
get_ni<- function(df, species, treat, boot_id){
  lambda <- with(df, lambda[ focal == species & treatment == treat & id == boot_id])[1]
  #print(lambda)
  gi <- with( df, g[focal == species & treatment == treat & id == boot_id])[1]
  #print(gi)
  si <- with( df, s[focal == species & treatment == treat & id == boot_id])[1]
  #print(si)
  
  ni<- ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}

get_ni(nls_boots, "ACWR", 1, "Bootstrap0002") #test

spp_treat_boot_combos$ni <- 0  # Add ni column to data frame

for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  treat <- spp_treat_boot_combos[i , "treatment"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  ni <- get_ni(nls_boots, sp1, treat, boot_id)
  spp_treat_boot_combos[i, "ni"] <- ni
  
}

#Get fitness differences ------

fitness_diff <- function(df1, df2, species1, species2, tr, boot_id) {
  
  ni <- with(df1, ni[focal == species1 & treatment == tr & id == boot_id])[1]
  #print(paste("ni: ",ni))
  nj <- with(df1, ni[focal == species2 & treatment == tr & id == boot_id])[1]
  #print(paste("nj: ",nj))
  
  aij <- with(df2, alpha[focal == species1 & competitor == species2 & treatment == tr & id == boot_id])
  #print(paste("aij: ",aij))
  
  aji <- with(df2, alpha[focal == species2 & competitor == species1 & treatment == tr & id == boot_id])
  #print(paste("aji: ",aji))
  
  ajj <- with(df2, alpha[focal == species2 & competitor == species2 & treatment == tr & id == boot_id])
  #print(paste("ajj: ",ajj))
  
  aii <- with(df2, alpha[focal == species1 & competitor == species1 & treatment == tr & id == boot_id])
  #print(paste("aii: ",aii))
  
  nn<- (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa<- sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij <- nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff(spp_treat_boot_combos, nls_boots, "HOMU", "URLI", 2, "Bootstrap0131") #test

with(spp_treat_boot_combos, ni[focal== "ACWR" & treatment == 1 & id == "Bootstrap0019"])

#add fitness difference column to data frame
spp_treat_boot_combos$fd <- 0

for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_boot_combos[i, "competitor"] %>% unlist
  trt <- spp_treat_boot_combos[i , "treatment"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  fd <- fitness_diff(spp_treat_boot_combos, nls_boots, sp1, sp2, trt, boot_id)
  spp_treat_boot_combos[i, "fd"] <- fd
  
}


#nls_boots <- read.csv("./output/nls_bootstraps_1000.csv")
#write_csv(spp_treat_boot_combos, "./output/spp_treat_boot_combos_1000_fullmodel.csv")
#spp_treat_boot_combos <- read.csv("./output/spp_treat_boot_combos_1000.csv")

##combine nls_boots (alpha and lambda data) with spp_treat_boot_combos (snd, ni, fd)

final_output_nls_boot <- merge(nls_boots, spp_treat_boot_combos, by = c("id", "focal", "competitor", "treatment"))
#write.csv(final_output_nls_boot, "./output/final_output_nls_boot_1000.csv")
#final_output_nls_boot <-read.csv("./output/final_output_nls_boot_1000.csv")
spp_treat_boot_combos <- final_output_nls_boot 


#Data frame with medians and sds for parameters----

spp_list <- sort(na.omit( unique(seed_data$focal)))
treat_list <- sort( na.omit( unique(seed_data$Tr)))
spp_treat_combos <- expand.grid(species = spp_list, treatment = treat_list)
comp_labels <- sort( na.omit( unique(seed_data$background) ))

spp_treat_boot_combos$sp_pair <- paste(spp_treat_boot_combos$focal, spp_treat_boot_combos$competitor, sep = "_")  

spp_treat_comp_combos <- expand.grid(focal = spp_list, competitor = comp_labels, treatment = treat_list)

nls_boot_pairs <- spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_sd = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0, 
         snd = 0, snd_low = 0, snd_high = 0, fd = 0, fd_low = 0, fd_high = 0)
nls_boot_pairs$sp_pair <- paste(nls_boot_pairs$focal, nls_boot_pairs$competitor, sep = "_")  


for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 <- spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt <- spp_treat_comp_combos[i , "treatment"] %>% unlist
  nls_boot_pairs[i, "alpha"] <- median(with(spp_treat_boot_combos, 
                                            alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                            ]), na.rm=TRUE)
  nls_boot_pairs[i, "alpha_sd"] <- sd(with(spp_treat_boot_combos, 
                                           alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                           ]), na.rm=TRUE)
  nls_boot_pairs[i, "alpha_low"] <- quantile(with(spp_treat_boot_combos, 
                                                  alpha[focal == sp1 & competitor == sp2 & 
                                                          treatment== treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "alpha_high"] <- quantile(with(spp_treat_boot_combos, 
                                                   alpha[focal == sp1 & competitor == sp2 & 
                                                           treatment== treatt]), 0.84, na.rm=TRUE)
  nls_boot_pairs[i, "lambda"] <- median(with(spp_treat_boot_combos, 
                                             lambda[focal == sp1 & competitor == sp2 & 
                                                      treatment == treatt]), na.rm=TRUE)
  nls_boot_pairs[i, "lambda_low"] <- quantile(with(spp_treat_boot_combos, 
                                                   lambda[focal == sp1 & competitor == sp2 & 
                                                            treatment == treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "lambda_high"] <- quantile(with(spp_treat_boot_combos, 
                                                    lambda[focal == sp1 & competitor == sp2 & 
                                                             treatment == treatt]), 0.84, na.rm=TRUE)
  nls_boot_pairs[i, "snd"] <- median(with(spp_treat_boot_combos, 
                                          snd[focal == sp1 & competitor == sp2 & 
                                                treatment == treatt]), na.rm=TRUE)
  nls_boot_pairs[i, "snd_low"] <- quantile(with(spp_treat_boot_combos, 
                                                snd[focal == sp1 & competitor == sp2 & 
                                                      treatment == treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "snd_high"] <- quantile(with(spp_treat_boot_combos, 
                                                 snd[focal == sp1 & competitor == sp2 & 
                                                       treatment == treatt]), 0.84, na.rm=TRUE)
  nls_boot_pairs[i, "fd"] <- median(with(spp_treat_boot_combos, 
                                         fd[focal == sp1 & competitor == sp2 & 
                                              treatment == treatt]),na.rm=TRUE)
  nls_boot_pairs[i, "fd_low"] <- quantile(with(spp_treat_boot_combos, 
                                               fd[focal == sp1 & competitor == sp2 & 
                                                    treatment == treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "fd_high"] <- quantile(with(spp_treat_boot_combos, 
                                                fd[focal == sp1 & competitor == sp2 & 
                                                     treatment == treatt]), 0.84, na.rm=TRUE)
}  

#write.csv(nls_boot_pairs, "./output/nls_boot_pairs_1000_full_model.csv")


