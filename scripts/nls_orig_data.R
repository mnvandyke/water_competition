### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

###This script uses non-linear least squares method to fit the annual plant model with the original collected data - No bootstrap
### and then calculates stabilizing niche and fitness differences for each pair under each treatment

library(tidyverse)
library(nlstools)
library(dplyr)
library(broom)
options(scipen = 5)

seed_data <- read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr <- ifelse(seed_data$treat == "W", 1, 2)
seed_data <- seed_data %>% replace_na(list(background = "ACWR")) # putting acwr instead of NA for the lambda individuals
#Need to get number of each species in the background for the model
seed_data$N_acwr <- ifelse(seed_data$background == "ACWR", seed_data$num_comp, 0)
seed_data$N_femi <- ifelse(seed_data$background == "FEMI", seed_data$num_comp, 0)
seed_data$N_homu <- ifelse(seed_data$background == "HOMU", seed_data$num_comp, 0)
seed_data$N_pler <- ifelse(seed_data$background == "PLER", seed_data$num_comp, 0)
seed_data$N_saco <- ifelse(seed_data$background == "SACO", seed_data$num_comp, 0)
seed_data$N_urli <- ifelse(seed_data$background == "URLI", seed_data$num_comp, 0)

acwr <- seed_data %>%
  filter(focal == "ACWR")


fit_test<- nls(log(num_seeds)~log(lambda[Tr]/(1+a_ACWR[Tr]*N_acwr+a_FEMI[Tr]*N_femi+a_HOMU[Tr]*N_homu+a_PLER[Tr]*N_pler+a_SACO[Tr]*N_saco+a_URLI[Tr]*N_urli)),
               data=acwr, start=list('lambda'= c(100,100), a_ACWR=c(0.1, 0.1), a_FEMI=c(0.1, 0.1), a_HOMU=c(0.1, 0.1), a_PLER=c(0.1, 0.1), a_SACO=c(0.1, 0.1),a_URLI=c(0.1, 0.1)),
               lower = c( 1, 1, rep(.001, 12) ), # setting alpha to be greater than .001  
               upper = c(10000, 10000, rep(2, 12)), 
               control = list(maxiter = 100000),
               algorithm = 'port') 
summary(fit_test)
estimates <- tidy(fit_test)

estimates$focal <- "ACWR"

estimates<-estimates%>% 
  select(term , estimate, focal) %>% 
  spread(term, estimate)

estimates <- estimates %>% pivot_longer(
  cols = starts_with("a"), 
  names_to = c("competitor","treatment"), 
  names_prefix = "a_", 
  names_sep = 4,
  values_to ="alpha"
)

estimates$lambda<-ifelse(estimates$treatment == 1, estimates$lambda1, estimates$lambda2)
estimates$lambda1<-NULL
estimates$lambda2<-NULL

#Now do that for all species and combine the data frames

spp_list <- sort(na.omit( unique(seed_data$focal)))

spp_combos <- expand.grid(species = spp_list)

# loop through the species list, fit the model, return the dataframe 
out <- list() #list for for loop output

for( i in 1:nrow( spp_combos )){ 
  temp_data <- seed_data %>%
    filter(focal == spp_combos[i,1])
  
  fit_test<- nls(log(num_seeds)~log(lambda[Tr]/(1+a_ACWR[Tr]*N_acwr+a_FEMI[Tr]*N_femi+a_HOMU[Tr]*N_homu+a_PLER[Tr]*N_pler+a_SACO[Tr]*N_saco+a_URLI[Tr]*N_urli)),
                 data=temp_data, start=list('lambda'= c(100,100), a_ACWR=c(0.1, 0.1), a_FEMI=c(0.1, 0.1), a_HOMU=c(0.1, 0.1), a_PLER=c(0.1, 0.1), a_SACO=c(0.1, 0.1),a_URLI=c(0.1, 0.1)),
                 lower = c( 1, 1, rep(.001, 12) ),# setting alpha to be greater than .001  
                 upper = c(10000, 10000, rep(2, 12)), 
                 control = list(maxiter = 100000),
                 algorithm = 'port') 
  
  df <- tidy(fit_test)
  
  df$focal <- spp_combos[i,1]
  
 df <- df %>% 
    select(term , estimate, focal) %>% 
    spread(term, estimate)
  
  df <- df %>% pivot_longer(
    cols = starts_with("a"), 
    names_to = c("competitor","treatment"), 
    names_prefix = "a_", 
    names_sep = 4,
    values_to ="alpha"
  )
  
  df$lambda<-ifelse(df$treatment == 1, df$lambda1, df$lambda2)
  df$lambda1<-NULL
  df$lambda2<-NULL
  
  out[[i]] <- df
}

all_fit <- do.call(rbind.data.frame, out)

##Calculate Stabilizing niche differences

stabilizing_niche_diff_func <- function(df, species1, species2, treat) {
  
  
  aij <- with(df, alpha[focal == species1 & competitor == species2 & treatment== treat])
  #print(paste("aij: ",aij))
  
  aji <- with(df, alpha[focal == species2 & competitor == species1 & treatment== treat])
  #print(paste("aji: ",aji))
  
  ajj <- with(df, alpha[focal == species2 & competitor == species2 & treatment== treat])
  #print(paste("ajj: ",ajj))
  
  aii <- with(df, alpha[focal == species1 & competitor == species1 & treatment== treat])
  #print(paste("aii: ",aii))
  
  snd <- (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

stabilizing_niche_diff_func(all_fit, "ACWR", "FEMI", 2) #test
all_fit$snd <- 0

for(i in 1:nrow(all_fit)) {
  
  sp1 <- all_fit[i, "focal"] %>% unlist
  sp2 <- all_fit[i, "competitor"] %>% unlist
  trt <- all_fit[i , "treatment"] %>% unlist
  snd <- stabilizing_niche_diff_func(all_fit, sp1, sp2, trt)
  all_fit[i, "snd"] <- snd
  
}

##Get ni-------
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 

s_g_data <- read.csv("./data/s_g_data.csv") #seed survival and germination data
all_fit <- merge(all_fit, s_g_data, by = "focal")
all_fit$X<- NULL


#ηi equation function

get_ni_func<- function(df, species, treat){
  lambda <- with(df, lambda[ focal == species & treatment == treat])[1]
  #print(lambda)
  gi <- with( df, g[focal == species & treatment == treat ])[1]
  #print(gi)
  si <- with( df, s[focal == species & treatment == treat])[1]
  #print(si)
  
  ni<- ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}


get_ni_func(all_fit, "ACWR", 1) #test
all_fit$ni <- 0  # Add ni column to data frame

for(i in 1:nrow(all_fit)) {
  
  sp1 <- all_fit[i, "focal"] %>% unlist
  trt <- all_fit[i , "treatment"] %>% unlist
  ni <- get_ni_func(all_fit, sp1, trt)
  all_fit[i, "ni"] <- ni
  
}

#Get fitness differences ------

fitness_diff_func <- function(df, species1, species2, treat) {
  
  ni <- with(df, ni[focal == species1 & treatment == treat])[1]
  #print(paste("ni: ",ni))
  nj <- with(df, ni[focal == species2 & treatment == treat])[1]
  #print(paste("nj: ",nj))
  
  aij <- with(df, alpha[focal == species1 & competitor == species2 & treatment == treat ])
  #print(paste("aij: ",aij))
  
  aji <- with(df, alpha[focal == species2 & competitor == species1 & treatment == treat])
  #print(paste("aji: ",aji))
  
  ajj <- with(df, alpha[focal == species2 & competitor == species2 & treatment == treat])
  #print(paste("ajj: ",ajj))
  
  aii <- with(df, alpha[focal == species1 & competitor == species1 & treatment == treat])
  #print(paste("aii: ",aii))
  
  nn<- (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa<- sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij <- nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff_func(all_fit, "FEMI", "HOMU", 1) #test

#add fitness difference column to data frame
all_fit$fd <- 0

for(i in 1:nrow(all_fit)) {
  
  sp1 <- all_fit[i, "focal"] %>% unlist
  sp2 <- all_fit[i, "competitor"] %>% unlist
  trt <- all_fit[i , "treatment"] %>% unlist
  fitdif <- fitness_diff_func(all_fit, sp1, sp2, trt)
  all_fit[i, "fd"] <- fitdif
}

all_fit$focal <- as.character(all_fit$focal)
all_fit$fd_superior <- ifelse(all_fit$fd < 1, 1/all_fit$fd, all_fit$fd) # identifying greater fitness difference
all_fit$fd_sup_sp <- ifelse(all_fit$fd <= 1, all_fit$focal, all_fit$competitor)

all_fit$coexist <- ifelse((all_fit$snd > (1-1/all_fit$fd_superior)), 1, 0 )
all_fit$sp_pair <- paste(all_fit$focal, all_fit$competitor, sep = "_")


