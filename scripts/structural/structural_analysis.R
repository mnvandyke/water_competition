### Analysis for "Small rainfall changes drive substantial changes in plant coexistence"
### Mary Van Dyke, mnvandyke@ucla.edu
### Last edit: 25 July 2022

### This file performs structural analysis for each species pair, triplet, quadruplet, quintuplet, and sextuplet.
### Adapted from code by Dr. Chris Johnson and from Saavedra et al. 2017 Ecological Monographs

library(dplyr)
options(scipen = 5)
source("./scripts/structural/toolbox_coexistence.R") # needs to run first to obtain needed functions
pars <- read.csv("./output/nls_boot_pairs_1000_full_model.csv")
#pars <- all_fit
sg <- read.csv("./data/s_g_data.csv")

pairs_struc <- pars[, c("focal", "competitor", "treatment", "alpha", "lambda") ]

structural_pair <- function(df, species1, species2, trt){
  
  pars_t <- df %>%
    filter(treatment == trt) 
  
sp1 <- species1
sp2 <- species2

if(sp1 == sp2) {return(c(0, 0, 0))} else
  
a11 <- pars_t$alpha[pars_t$focal == sp1 & pars_t$competitor == sp1]*sg$g[sg$focal== sp1]
a12 <- pars_t$alpha[pars_t$focal == sp1 & pars_t$competitor == sp2]*sg$g[sg$focal== sp2]

a22 <- pars_t$alpha[pars_t$focal == sp2 & pars_t$competitor == sp2]*sg$g[sg$focal== sp2]
a21 <- pars_t$alpha[pars_t$focal == sp2 & pars_t$competitor == sp1]*sg$g[sg$focal== sp1]

r1 <- ((sg$g[sg$focal==sp1]*pars_t$lambda[pars_t$focal == sp1][1])/(1 - (1-sg$g[sg$focal==sp1])*sg$s[sg$focal==sp1])) - 1
r2 <- ((sg$g[sg$focal==sp2]*pars_t$lambda[pars_t$focal == sp2][1])/(1 - (1-sg$g[sg$focal==sp2])*sg$s[sg$focal==sp2])) - 1

alphas <- matrix(c(a11,a12,a21,a22), nrow=2, ncol=2, byrow = TRUE)
r <- c(r1,r2)
om <- 10^Omega(alphas) # structural analog of the niche difference (on linear scale)
th <- theta(alphas, r) # structural analog of the fitness difference
comp<- test_feasibility(alphas,r) # test feasibility of entire community (1 = feasible)
test_feasibility_pairs(alphas,r) # test feasibility of each pair (1 = feasible)

return(c(om, th, comp))}


structural_pair(pairs_struc, "URLI", "FEMI", 2) #test

pairs_struc$omega <- 0
pairs_struc$theta <- 0
pairs_struc$outcome <- 0

for( i in 1:nrow(pairs_struc)) {
  focal1 <- pairs_struc[i, "focal"] %>% unlist
  competitor2 <- pairs_struc[i, "competitor"] %>% unlist
  treatt <- pairs_struc[i , "treatment"] %>% unlist
  
  struc <- structural_pair(pairs_struc, focal1, competitor2, treatt)
  print(struc)
  
  pairs_struc[i, "omega"] <- struc[1] %>% unlist
  pairs_struc[i, "theta"] <- struc[2] %>% unlist
  pairs_struc[i, "outcome"] <- struc[3] %>% unlist
  
}

pairs_struc$species <- paste(substr(pairs_struc$focal, start = 1, stop = 2), 
                            substr(pairs_struc$competitor, start = 1, stop = 2), 
                                   sep = "-")
pairs_struc$no_sp <- 2

pairs_struc$test <- pairs_struc$omega *45
                                                                

##Triplets

structural_trip <- function(df, species1, species2, species3, trt){
pars3 <- df %>%
  filter(treatment == trt) 
sp1 <- species1
sp2 <- species2
sp3 <- species3

if(sp1 == sp2) {return(c(0, 0, 0))} else
  if(sp1 == sp3) {return(c(0, 0, 0))} else
    if(sp2 == sp3) {return(c(0, 0, 0))} else
  
a11 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a12 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a13 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]

a22 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a21 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a23 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]

a31 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a32 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a33 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]

r1 <- ((sg$g[sg$focal==sp1]*pars3$lambda[pars3$focal == sp1][1])/(1 - (1-sg$g[sg$focal==sp1])*sg$s[sg$focal==sp1])) - 1
r2 <- ((sg$g[sg$focal==sp2]*pars3$lambda[pars3$focal == sp2][1])/(1 - (1-sg$g[sg$focal==sp2])*sg$s[sg$focal==sp2])) - 1
r3 <- ((sg$g[sg$focal==sp3]*pars3$lambda[pars3$focal == sp3][1])/(1 - (1-sg$g[sg$focal==sp3])*sg$s[sg$focal==sp3])) - 1


# Outputs
alpha <- matrix(c(a11,a12,a13,a21,a22,a23,a31,a32,a33), nrow=3, ncol=3, byrow = TRUE)
r <- c(r1,r2,r3)
om <- 10^Omega(alpha) # structural analog of the niche difference (on linear scale)
th<-theta(alpha,r) # structural analog of the fitness difference
out<-test_feasibility(alpha,r) # test feasibility of entire community (1 = feasible)
overlap <- compute_overlap(alpha,10000)

return(c(om, th, out, overlap))
}

structural_trip(pars, "FEMI", "URLI", "PLER", 2)

spp_list <- sort(na.omit( unique(pars$focal)))
pars_trip <- pars[, c("focal", "competitor", "treatment", "alpha", "lambda") ]


triplets<-as.matrix(t(combn(unique(spp_list),3)))
triplets <-data.frame(triplets[,1], triplets[,2], triplets[,3])
triplets <- merge(triplets, c(1,2))
colnames(triplets) <- c("species_1", "species_2", "species_3", "treatment")
triplets$omega <- 0
triplets$theta <- 0
triplets$outcome <- 0



for( i in 1:nrow(triplets)) {
  sp_1 <- triplets[i, "species_1"] %>% unlist
  sp_2 <- triplets[i, "species_2"] %>% unlist
  sp_3 <- triplets[i , "species_3"] %>% unlist
  treatt <- triplets[i, "treatment"] %>% unlist
  
  struc <- structural_trip(pairs_struc, sp_1, sp_2, sp_3, treatt)
  print(struc)
  
  triplets[i, "omega"] <- struc[1] %>% unlist
  triplets[i, "theta"] <- struc[2] %>% unlist
  triplets[i, "outcome"] <- struc[3] %>% unlist
  
}

triplets$species <- paste(substr(triplets$species_1, start = 1, stop = 2), substr(triplets$species_2, start = 1, stop = 2),substr(triplets$species_3, start = 1, stop = 2), sep="-")
triplets$no_sp <- 3


# QUADRUPLETS (use parameters for triplets above and add 4th species)-----
# parameters

structural_4 <- function(df, species1, species2, species3,species4, trt){
  pars3 <- df %>%
    filter(treatment == trt) 
  sp1 <- species1
  sp2 <- species2
  sp3 <- species3
  sp4 <- species4
  
  if(sp1 == sp2) {return(c(0, 0, 0))} else
    if(sp1 == sp3) {return(c(0, 0, 0))} else
      if(sp1 == sp4) {return(c(0, 0, 0))} else
        if(sp2 == sp3) {return(c(0, 0, 0))} else
          if(sp2 == sp4) {return(c(0, 0, 0))} else
            if(sp3 == sp4) {return(c(0, 0, 0))} else

a11 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a12 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a13 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a14 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]

a22 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a21 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a23 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a24 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]

a31 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a32 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a33 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a34 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]

a41 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a42 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a43 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a44 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]

r1 <- ((sg$g[sg$focal==sp1]*pars3$lambda[pars3$focal == sp1][1])/(1 - (1-sg$g[sg$focal==sp1])*sg$s[sg$focal==sp1])) - 1
r2 <- ((sg$g[sg$focal==sp2]*pars3$lambda[pars3$focal == sp2][1])/(1 - (1-sg$g[sg$focal==sp2])*sg$s[sg$focal==sp2])) - 1
r3 <- ((sg$g[sg$focal==sp3]*pars3$lambda[pars3$focal == sp3][1])/(1 - (1-sg$g[sg$focal==sp3])*sg$s[sg$focal==sp3])) - 1
r4 <- ((sg$g[sg$focal==sp4]*pars3$lambda[pars3$focal == sp4][1])/(1 - (1-sg$g[sg$focal==sp4])*sg$s[sg$focal==sp4])) - 1

# Outputs
alpha <- matrix(c(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44), nrow=4, ncol=4, byrow = TRUE)
r <- c(r1,r2,r3,r4)
om <- 10^Omega(alpha) # structural analog of the niche difference (on linear scale)
th <- theta(alpha,r) # structural analog of the fitness difference
outcome <- test_feasibility(alpha,r) # test feasibility of entire community (1 = feasible)

return(c(om, th, outcome))
}

structural_4(pairs_struc, "ACWR", "FEMI", "HOMU", "PLER", 1)

quarts <- as.matrix(t((combn(unique(spp_list),4))))
quarts <- data.frame(quarts[,1], quarts[,2], quarts[,3], quarts[,4])
quarts <- merge(quarts, c(1,2))
colnames(quarts) <- c("species_1", "species_2", "species_3", "species_4", "treatment")
quarts$omega <- 0
quarts$theta <- 0
quarts$outcome <- 0



for( i in 1:nrow(quarts)) {
  sp_1 <- quarts[i, "species_1"] %>% unlist
  sp_2 <- quarts[i, "species_2"] %>% unlist
  sp_3 <- quarts[i , "species_3"] %>% unlist
  sp_4 <- quarts[i , "species_4"] %>% unlist
  treatt <- quarts[i, "treatment"] %>% unlist
  
  struc <- structural_4(pairs_struc, sp_1, sp_2, sp_3, sp_4, treatt)
  print(struc)
  
  quarts[i, "omega"] <- struc[1] %>% unlist
  quarts[i, "theta"] <- struc[2] %>% unlist
  quarts[i, "outcome"] <- struc[3] %>% unlist
  
}

quarts$species <- quarts$species <- paste(substr(quarts$species_1, start = 1, stop = 2), 
                                          substr(quarts$species_2, start = 1, stop = 2),
                                          substr(quarts$species_3, start = 1, stop = 2), 
                                          substr(quarts$species_4, start = 1, stop = 2),
                                          sep="-")


quarts$no_sp <- 4

##FIVES

structural_5 <- function(df, species1, species2, species3,species4, species5, trt){
  pars3 <- df %>%
    filter(treatment == trt) 
  sp1 <- species1
  sp2 <- species2
  sp3 <- species3
  sp4 <- species4
  sp5 <- species5
  
  if(sp1 == sp2) {return(c(0, 0, 0))} else
    if(sp1 == sp3) {return(c(0, 0, 0))} else
      if(sp1 == sp4) {return(c(0, 0, 0))} else
        if(sp1 == sp5) {return(c(0, 0, 0))} else
          if(sp2 == sp3) {return(c(0, 0, 0))} else
            if(sp2 == sp4) {return(c(0, 0, 0))} else
              if(sp2 == sp5) {return(c(0, 0, 0))} else
                if(sp3 == sp4) {return(c(0, 0, 0))} else
                  if(sp4 == sp5) {return(c(0, 0, 0))} else
              
  a11 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
  a12 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
  a13 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
  a14 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
  a15 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
  
  a22 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
  a21 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
  a23 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
  a24 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
  a25 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
  
  
  a31 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
  a32 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
  a33 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
  a34 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
  a35 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
  
  
  a41 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
  a42 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
  a43 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
  a44 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
  a45 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
  
  a51 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
  a52 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
  a53 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
  a54 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
  a55 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
  
  
  r1 <- ((sg$g[sg$focal==sp1]*pars3$lambda[pars3$focal == sp1][1])/(1 - (1-sg$g[sg$focal==sp1])*sg$s[sg$focal==sp1])) - 1
  r2 <- ((sg$g[sg$focal==sp2]*pars3$lambda[pars3$focal == sp2][1])/(1 - (1-sg$g[sg$focal==sp2])*sg$s[sg$focal==sp2])) - 1
  r3 <- ((sg$g[sg$focal==sp3]*pars3$lambda[pars3$focal == sp3][1])/(1 - (1-sg$g[sg$focal==sp3])*sg$s[sg$focal==sp3])) - 1
  r4 <- ((sg$g[sg$focal==sp4]*pars3$lambda[pars3$focal == sp4][1])/(1 - (1-sg$g[sg$focal==sp4])*sg$s[sg$focal==sp4])) - 1
  r5 <- ((sg$g[sg$focal==sp5]*pars3$lambda[pars3$focal == sp5][1])/(1 - (1-sg$g[sg$focal==sp5])*sg$s[sg$focal==sp5])) - 1
  
  # Outputs
  alphas <- matrix(c(a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55), nrow=5, ncol=5, byrow = TRUE)
  r <- c(r1,r2,r3,r4,r5)
  om <- 10^Omega(alphas) # structural analog of the niche difference (on linear scale)
  th <- theta(alphas,r) # structural analog of the fitness difference
  outcome <- test_feasibility(alphas,r) # test feasibility of entire community (1 = feasible)
  
  return(c(om, th, outcome))
}

structural_5(pairs_struc, "ACWR", "FEMI", "HOMU", "PLER", "SACO", 1)

quints <- as.matrix(t((combn(unique(spp_list),5))))
quints <- data.frame(quints[,1], quints[,2], quints[,3], quints[,4], quints[,5])
quints <- merge(quints, c(1,2))
colnames(quints) <- c("species_1", "species_2", "species_3", "species_4", "species_5", "treatment")
quints$omega <- 0
quints$theta <- 0
quints$outcome <- 0



for( i in 1:nrow(quints)) {
  sp_1 <- quints[i, "species_1"] %>% unlist
  sp_2 <- quints[i, "species_2"] %>% unlist
  sp_3 <- quints[i , "species_3"] %>% unlist
  sp_4 <- quints[i , "species_4"] %>% unlist
  sp_5 <- quints[i , "species_5"] %>% unlist
  treatt <- quints[i, "treatment"] %>% unlist
  
  struc <- structural_5(pairs_struc, sp_1, sp_2, sp_3, sp_4, sp_5, treatt)
  print(struc)
  
  quints[i, "omega"] <- struc[1] %>% unlist
  quints[i, "theta"] <- struc[2] %>% unlist
  quints[i, "outcome"] <- struc[3] %>% unlist
  
}

quints$species <- quints$species <- paste(substr(quints$species_1, start = 1, stop = 2), 
                                          substr(quints$species_2, start = 1, stop = 2),
                                          substr(quints$species_3, start = 1, stop = 2), 
                                          substr(quints$species_4, start = 1, stop = 2),
                                          substr(quints$species_5, start = 1, stop = 2),
                                          sep="-")
quints$no_sp <- 5

##Sixtuplet

structural_6 <- function(df, species1, species2, species3,species4, species5, species6, trt){
  pars3 <- df %>%
    filter(treatment == trt) 
  
  sp1 <- species1
  sp2 <- species2
  sp3 <- species3
  sp4 <- species4
  sp5 <- species5
  sp6 <- species6
  
  if(sp1 == sp2) {return(c(0, 0, 0))} else
    if(sp1 == sp3) {return(c(0, 0, 0))} else
      if(sp1 == sp4) {return(c(0, 0, 0))} else
        if(sp1 == sp5) {return(c(0, 0, 0))} else
          if(sp2 == sp3) {return(c(0, 0, 0))} else
            if(sp2 == sp4) {return(c(0, 0, 0))} else
              if(sp2 == sp5) {return(c(0, 0, 0))} else
                if(sp3 == sp4) {return(c(0, 0, 0))} else
                  if(sp4 == sp5) {return(c(0, 0, 0))} else
                    if(sp1 == sp6) {return(c(0, 0, 0))} else
                      if(sp2 == sp6) {return(c(0, 0, 0))} else
                        if(sp3 == sp6) {return(c(0, 0, 0))} else
                          if(sp4 == sp6) {return(c(0, 0, 0))} else
                            if(sp5 == sp6) {return(c(0, 0, 0))} else

a11 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a12 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a13 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a14 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a15 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a16 <- pars3$alpha[pars3$focal == sp1 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

a22 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a21 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a23 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a24 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a25 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a26 <- pars3$alpha[pars3$focal == sp2 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

a31 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a32 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a33 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a34 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a35 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a36 <- pars3$alpha[pars3$focal == sp3 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

a41 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a42 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a43 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a44 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a45 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a46 <- pars3$alpha[pars3$focal == sp4 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

a51 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a52 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a53 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a54 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a55 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a56 <- pars3$alpha[pars3$focal == sp5 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

a61 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp1]*sg$g[sg$focal== sp1]
a62 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp2]*sg$g[sg$focal== sp2]
a63 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp3]*sg$g[sg$focal== sp3]
a64 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp4]*sg$g[sg$focal== sp4]
a65 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp5]*sg$g[sg$focal== sp5]
a66 <- pars3$alpha[pars3$focal == sp6 & pars3$competitor == sp6]*sg$g[sg$focal== sp6]

r1 <- ((sg$g[sg$focal==sp1]*pars3$lambda[pars3$focal == sp1][1])/(1 - (1-sg$g[sg$focal==sp1])*sg$s[sg$focal==sp1])) - 1
r2 <- ((sg$g[sg$focal==sp2]*pars3$lambda[pars3$focal == sp2][1])/(1 - (1-sg$g[sg$focal==sp2])*sg$s[sg$focal==sp2])) - 1
r3 <- ((sg$g[sg$focal==sp3]*pars3$lambda[pars3$focal == sp3][1])/(1 - (1-sg$g[sg$focal==sp3])*sg$s[sg$focal==sp3])) - 1
r4 <- ((sg$g[sg$focal==sp4]*pars3$lambda[pars3$focal == sp4][1])/(1 - (1-sg$g[sg$focal==sp4])*sg$s[sg$focal==sp4])) - 1
r5 <- ((sg$g[sg$focal==sp5]*pars3$lambda[pars3$focal == sp5][1])/(1 - (1-sg$g[sg$focal==sp5])*sg$s[sg$focal==sp5])) - 1
r6 <- ((sg$g[sg$focal==sp6]*pars3$lambda[pars3$focal == sp6][1])/(1 - (1-sg$g[sg$focal==sp6])*sg$s[sg$focal==sp6])) - 1

# Outputs
alphas <- matrix(c(a11,a12,a13,a14,a15,a16, a21,a22,a23,a24,a25,a26, a31,a32,a33,a34,a35,a36, a41,a42,a43,a44,a45,a46, a51,a52,a53,a54,a55,a56, a61,a62,a63,a64,a65,a66), nrow=6, ncol=6, byrow = TRUE)
r <- c(r1,r2,r3,r4, r5,r6)
om<-10^Omega(alphas) # structural analog of the niche difference (on linear scale)
th <-theta(alphas,r) # structural analog of the fitness difference
outcome<- test_feasibility(alphas,r) # test feasibility of entire community (1 = feasible)

return(c(om, th, outcome))
}

structural_6(pairs_struc, "ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI", 1)

sexts <- as.matrix(t((combn(unique(spp_list),6))))
sexts <- data.frame(sexts[,1], sexts[,2], sexts[,3], sexts[,4], sexts[,5], sexts[,6])
sexts <- merge(sexts, c(1,2))
colnames(sexts) <- c("species_1", "species_2", "species_3", "species_4", "species_5", "species_6", "treatment")
sexts$omega <- 0
sexts$theta <- 0
sexts$outcome <- 0

for( i in 1:nrow(sexts)) {
  sp_1 <- sexts[i, "species_1"] %>% unlist
  sp_2 <- sexts[i, "species_2"] %>% unlist
  sp_3 <- sexts[i , "species_3"] %>% unlist
  sp_4 <- sexts[i , "species_4"] %>% unlist
  sp_5 <- sexts[i , "species_5"] %>% unlist
  sp_6 <- sexts[i , "species_6"] %>% unlist
  treatt <- sexts[i, "treatment"] %>% unlist
  
  struc <- structural_6(pairs_struc, sp_1, sp_2, sp_3, sp_4, sp_5, sp_6, treatt)
  print(struc)
  
  sexts[i, "omega"] <- struc[1] %>% unlist
  sexts[i, "theta"] <- struc[2] %>% unlist
  sexts[i, "outcome"] <- struc[3] %>% unlist
  
}

sexts$species <- sexts$species <- paste(substr(sexts$species_1, start = 1, stop = 2), 
                                          substr(sexts$species_2, start = 1, stop = 2),
                                          substr(sexts$species_3, start = 1, stop = 2), 
                                          substr(sexts$species_4, start = 1, stop = 2),
                                          substr(sexts$species_5, start = 1, stop = 2),
                                          substr(sexts$species_6, start = 1, stop = 2),
                                          sep="-")
sexts$no_sp <- 6


##Output Combined Data frame

#put all data frames into list
pairs_struc <- pairs_struc %>%
  filter(focal != competitor)
df_list <- list(pairs_struc, triplets, quarts, quints, sexts)      

#merge all data frames together
merged<-df_list %>% reduce(full_join, by=c("species", "treatment", "omega", "theta", "outcome", "no_sp"))
merged <- merged[,c("species", "treatment", "omega", "theta", "outcome", "no_sp")]

write.csv(merged, "./output/structural_analysis_output.csv")

##Simulation of Triplets to check coexistence with lambdas and alphas\][]
s_g_data <- read.csv("./data/s_g_data.csv")

getNs <- function(N1, N2, N3){
  N0i <- N1
  N0j <- N2
  N0k <- N3

lambda_i = with(pars, lambda[focal == sp_i & treatment== trt])[1]
lambda_j = with(pars, lambda[focal == sp_j & treatment== trt])[1]
lambda_k = with(pars, lambda[focal == sp_k & treatment== trt])[1]
aii = with(pars, alpha[focal == sp_i & competitor == sp_i & treatment== trt])[1]
aij = with(pars, alpha[focal == sp_i & competitor == sp_j & treatment== trt])[1]
aik = with(pars, alpha[focal == sp_i & competitor == sp_k & treatment== trt])[1]
ajj = with(pars, alpha[focal == sp_j & competitor == sp_j & treatment== trt])[1]
aji = with(pars, alpha[focal == sp_j & competitor == sp_i & treatment== trt])[1]
ajk = with(pars, alpha[focal == sp_j & competitor == sp_k & treatment== trt])[1]
akk = with(pars, alpha[focal == sp_k & competitor == sp_k & treatment== trt])[1]
aki = with(pars, alpha[focal == sp_k & competitor == sp_i & treatment== trt])[1]
akj = with(pars, alpha[focal == sp_k & competitor == sp_j & treatment== trt])[1]
gi = with(s_g_data, g[focal == sp_i])[1]
gj = with(s_g_data, g[focal == sp_j])[1]
gk = with(s_g_data, g[focal == sp_k])[1]
si = with(s_g_data, s[focal == sp_i])[1]
sj = with(s_g_data, s[focal == sp_j])[1]
sk = with(s_g_data, s[focal == sp_k])[1]


Ni <- ((1-gi)*si + gi*(lambda_i/(1 + aii*gi*N0i + aij*gj*N0j + aik*gk*N0k)))*N0i
Nj <- ((1-gj)*sj + gj*(lambda_j/(1 + ajj*gj*N0j + aji*gi*N0i + ajk*gk*N0k)))*N0j 
Nk <- ((1-gk)*sk + gk*(lambda_k/(1 + akk*gk*N0k + akj*gj*N0j + aki*gi*N0i)))*N0k

return(c(Ni, Nj, Nk))

}

sp_i = "SACO"
sp_j = "URLI"
sp_k = "PLER"
trt = 2

getNs(100, 100, 100)

sim <- data.frame(time=1:1000,
                  Ni = 0,
                  Nj = 0,
                  Nk = 0)
  
sim[1, "Ni"]<- 10
sim[1, "Nj"]<- 1000
sim[1, "Nk"]<- 1000

for(t in 1:nrow(sim)){
  Nns <- getNs(sim[t, "Ni"], sim[t, "Nj"], sim[t, "Nk"])

  sim[(t+1), "Ni"]<- Nns[1]
  sim[(t+1), "Nj"]<- Nns[2]
  sim[(t+1), "Nk"]<- Nns[3]
  
}

ggplot(sim)+
  geom_point(aes(time, Ni), color = "red")+
  geom_point(aes(time, Nj), color = "blue")+
  geom_point(aes(time, Nk), color = "green")


