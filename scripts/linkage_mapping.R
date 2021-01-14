##### INTRODUCTION TO LINKAGE MAPPING #####

library(AlphaSimR)
library(tidyverse)

# 1. INTRODUCTION TO LINKAGE MAPPING #####

n_chr <- 1 # chromosomes
n_seg <- 10 # segregation sites
n_off <- 100 # number of offspring per cross

#___1.1 Develop genetic map and define haplotypes ####
genMap <- list(seq(0, 1, length.out = n_seg)) # list of segregation sites
genMap

chr1 <- c(rep(1, n_seg), # haplotype 1 of heterozygote parent
          rep(0, n_seg), # haplotype 2 of heterozygote parent
          rep(0, n_seg), # haplotype 1 of homozygote parent
          rep(0, n_seg)) # haplotype 2 of homozygote parent)

chr1 <- matrix(chr1, nrow = 4, ncol = n_seg, byrow = TRUE) # coerce into a matrix
chr1

pop_haplos <- newMapPop(genMap = genMap,
                        haplotypes = list(chr1))

#___1.2 Set simulation parameters and create parental population ####

SP = SimParam$new(pop_haplos) # create a variable for new simulation parameters
SP$setSexes("yes_sys") # set sexes systematically

pop_0 <- newPop(rawPop = pop_haplos, # object of MapPop-class
                simParam = SP) # SimParam objects

#___1.3 Cross parents and create offspring population ####

crossPlan <- matrix(c(1, 2),
                    nrow = 1,
                    ncol = 2) # develop a simple cross plan
crossPlan

pop_1 <- makeCross(pop = pop_0, # Pop-class object
                   crossPlan = crossPlan, # two column matrix indicating female/male cross made
                   nProgeny = n_off, # number of progeny per cross
                   simParam = SP) # simParam object

#___1.4 Pull genotypes of parents and offspring ####

df_P <- pullSegSiteGeno(pop = pop_0, # Pop-class object
                      chr = 1, # chromosome number, NULL = all
                      simParam = SP) # SimParam object
df_P <- as.data.frame(df_P)

df_O <- pullSegSiteGeno(pop = pop_1, # Pop-class object
                          chr = 1, # chromosome number, NULL = all
                          simParam = SP) # SimParam object
df_O <- as.data.frame(df_O)
df_O

#___1.5 Pairwise analysis for recombinant genotype frequency #####

# develop a pairwise plan

a <- 1:(n_seg-1)
ss1 <- rep(a, (n_seg-a))
foo <- function(b){c(1:(n_seg))[b:n_seg]}
ss2 <- unlist(lapply(2:(n_seg), FUN = foo))
pair_plan <- data.frame(ss1, ss2)  
pair_plan      

# develop a function for calculating recombinants based on pairs

pws <- 1:dim(pair_plan)[1]

rgf <- function(pws){
  
  ssA <- pair_plan[pws, 1]
  ssB <- pair_plan[pws, 2]
  
  p1 <- paste(df_P[1,ssA], df_P[2,ssB], sep = "")
  p2 <- paste(df_P[2,ssA], df_P[1,ssB], sep = "")
  
  dplyr::select(df_O, df_O[,ssA], df_O[,ssB]) %>%
    mutate(pair = paste(df_O[,ssA], df_O[,ssB], sep = "")) %>%
    mutate(recomb = ifelse(pair %in% c(p1, p2), 1, 0)) %>%
    summarize(sum = sum(recomb)/length(pop_1@id))  
}

reco <- unlist(lapply(pws, FUN = rgf))
pair_plan <- pair_plan %>%
  mutate(reco = reco)

# actual linkage distances between segregation sites
pair_plan %>%
  group_by(ss1) %>%
  mutate(min = min(reco)) %>%
  ungroup() %>%
  filter(reco == min)

# verify the order of the segregation sites
pair_plan %>%
  filter(ss1 == 1)





# 2. LINKAGE MAPPING USING THE ASMAP PACKAGE ####