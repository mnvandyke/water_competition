rm(list = ls())

library(tidyverse)
library(nlstools)
library(dplyr)
library(rsample)
library(broom)
library(purrr)
library(ggplot2)
library(ggrepel)
options(scipen = 5)
set.seed(3)

m = 6 # m is the number of competitors 
seed_data <- read.csv("./data/drought_seed_production_data.csv")

# function to find lambda and alpha parameter estimates, the split is a component of the bootstrap 
fit_nls_on_bootstrap <- function(split) {
  tryCatch({
    nls( y ~ log10(lambda/(1 + N*a[comp] )), 
         analysis(split), 
         start = list('lambda' = 200, 'a' = rep(0.1, m)), # these are initial values
         lower = c( 1, rep(.001, m) ),# forcing alpha to be greater than .001  
         #.002 is 2.5% of all alphas when lower bound was 0.001 but having alphas that small made 
         # some URLI stabilizing niche differences unrealistically large >>>1
         upper = c(10000, rep(2, m)), 
         control = list(maxiter = 100000),
         algorithm = 'port')  }, 
    error = function(e) { #Make it keep trying if it fails to converge at first
      print("failed to converge")  
    })
}
# Try different starting values ex: lambda = mean(seed_data$num_seeds)

#Create a data frame with all combinations of treatment and species 
spp_list <- sort(na.omit( unique(seed_data$focal)))

treat_list <- sort( na.omit( unique(seed_data$treat)))
spp_treat_combos <- expand.grid(species = spp_list, treat = treat_list)

# loop through the species list, fit the model, return the dataframe 
out <- list() #list for for loop output

for( i in 1:nrow( spp_treat_combos )){ 
  temp_data <- seed_data %>%
    filter(focal == spp_treat_combos[i,1]) %>%
    filter(treat == spp_treat_combos[i,2]) %>%
    mutate(N = num_comp, y = log10(num_seeds))
  temp_data$comp <- factor(temp_data$background) 
  temp_data$comp[is.na(temp_data$comp)] <- "PLER" #This is putting a competitor name in the lambda plots
  #b/c it's NA, Species identity shouldnt change anything
  
  #print((temp_data)) 
  
  boots <- bootstraps(temp_data, times = 1000)
  
  
  boot_models <- boots %>%
    mutate(model = map(splits, fit_nls_on_bootstrap),
           coef_info = map(model, tidy))
  
  boot_coefs <- boot_models %>%
    unnest(coef_info)
  
  boot_coefs$focal <- spp_treat_combos[i, 1]
  boot_coefs$treatment <- spp_treat_combos[i, 2]
  #print(head(boot_coefs))
  
  #This is removing a random column that sometimes appears
  out[[i]] <- boot_coefs
  if(!is.null(out[[i]][["x"]])) {
    out[[i]][["x"]]<- NULL
  }
}


test <- do.call(rbind.data.frame, out)

comp_labels <- sort( na.omit( unique(seed_data$background) ))
comp_labels
nls_boots<- test %>% 
  select( id, term , estimate, focal, treatment) %>% 
  spread(term, estimate) %>% 
  gather( competitor, alpha, starts_with('a')) %>% 
  mutate( competitor = factor( competitor, labels = comp_labels))  %>% 
  group_by( id)
#write_csv(nls_boots, "./output/nls_bootstraps_1000.csv")
#nls_boots<-read.csv("./output/nls_bootstraps_1000.csv")

#How many alphas are 0.001?
nrow(nls_boots[nls_boots$alpha == 0.001,])/72000

boot_list <- sort( na.omit( unique(nls_boots$id)))
spp_treat_boot_combos <- expand.grid(focal = spp_list, treat = treat_list, competitor = comp_labels, id = boot_list)

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
stabilizing_niche_diff(nls_boots, "ACWR", "URLI", "D", "Bootstrap0005")

#Add snd column to data frame
spp_treat_boot_combos$snd <- 0
spp_treat_boot_combos$treatment <- spp_treat_boot_combos$treat

for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_boot_combos[i, "competitor"] %>% unlist
  treatment <- spp_treat_boot_combos[i , "treat"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  snd <- stabilizing_niche_diff(nls_boots, sp1, sp2, treatment, boot_id)
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


get_ni(nls_boots, "ACWR", "D", "Bootstrap0001") #test
spp_treat_boot_combos$ni <- 0  # Add ni column to data frame

for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  treatment <- spp_treat_boot_combos[i , "treat"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  ni <- get_ni(nls_boots, sp1, treatment, boot_id)
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

#test
fitness_diff(spp_treat_boot_combos, nls_boots, "HOMU", "URLI", "W", "Bootstrap0131") #test

with(spp_treat_boot_combos, ni[focal== "ACWR" & treatment == "W" & id == "Bootstrap0019"])
#add fitness difference column to data frame
spp_treat_boot_combos$fd <- 0

for(i in 1:nrow(spp_treat_boot_combos)) {
  
  sp1 <- spp_treat_boot_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_boot_combos[i, "competitor"] %>% unlist
  trt <- spp_treat_boot_combos[i , "treat"] %>% unlist
  boot_id <- spp_treat_boot_combos[i, "id"] %>% unlist
  fd <- fitness_diff(spp_treat_boot_combos, nls_boots, sp1, sp2, trt, boot_id)
  spp_treat_boot_combos[i, "fd"] <- fd
  
}


#nls_boots <- read.csv("./output/nls_bootstraps_1000.csv")
#write_csv(spp_treat_boot_combos, "./output/spp_treat_boot_combos_1000.csv")
#spp_treat_boot_combos <- read.csv("./output/spp_treat_boot_combos_1000.csv")

##combine nls_boots (alpha and lambda data) with spp_treat_boot_combos (snd, ni, fd)

final_output_nls_boot <- merge(nls_boots, spp_treat_boot_combos, by = c("id", "focal", "competitor", "treatment"))
#write.csv(final_output_nls_boot, "./output/final_output_nls_boot_1000.csv")

spp_treat_boot_combos <- final_output_nls_boot 


#Data frame with medians and sds for parameters----

spp_list <- sort(na.omit( unique(seed_data$focal)))
treat_list <- sort( na.omit( unique(seed_data$treat)))
spp_treat_combos <- expand.grid(species = spp_list, treat = treat_list)
comp_labels <- sort( na.omit( unique(seed_data$background) ))

spp_treat_boot_combos$sp_pair <- paste(spp_treat_boot_combos$focal, spp_treat_boot_combos$competitor, sep = "_")  

spp_treat_comp_combos <- expand.grid(focal = spp_list,competitor = comp_labels, treat = treat_list)

nls_boot_pairs <- spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0, 
         snd = 0, snd_low = 0, snd_high = 0, fd = 0, fd_low = 0, fd_high = 0)
nls_boot_pairs$sp_pair <- paste(nls_boot_pairs$focal, nls_boot_pairs$competitor, sep = "_")  


for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 <- spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 <- spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt <- spp_treat_comp_combos[i , "treat"] %>% unlist
  nls_boot_pairs[i, "alpha"] <- median(with(spp_treat_boot_combos, 
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
                                                treat == treatt]), na.rm=TRUE)
  nls_boot_pairs[i, "snd_low"] <- quantile(with(spp_treat_boot_combos, 
                                                snd[focal == sp1 & competitor == sp2 & 
                                                      treat == treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "snd_high"] <- quantile(with(spp_treat_boot_combos, 
                                                 snd[focal == sp1 & competitor == sp2 & 
                                                       treat == treatt]), 0.84, na.rm=TRUE)
  nls_boot_pairs[i, "fd"] <- median(with(spp_treat_boot_combos, 
                                         fd[focal == sp1 & competitor == sp2 & 
                                              treat == treatt]),na.rm=TRUE)
  nls_boot_pairs[i, "fd_low"] <- quantile(with(spp_treat_boot_combos, 
                                               fd[focal == sp1 & competitor == sp2 & 
                                                    treat == treatt]), 0.16, na.rm=TRUE)
  nls_boot_pairs[i, "fd_high"] <- quantile(with(spp_treat_boot_combos, 
                                                fd[focal == sp1 & competitor == sp2 & 
                                                     treat == treatt]), 0.84, na.rm=TRUE)
}  

#write.csv(nls_boot_pairs, "./output/nls_boot_pairs_1000.csv")
