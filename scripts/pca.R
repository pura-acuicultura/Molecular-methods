##### DIFFERENTIATION OF POPULATIONS USING PCA #####

library(AlphaSimR)
library(tidyverse)

# 1. CREATE FOUNDER HAPLOTYPES #####
#__ 1.1 Founder haplotypes with quickHaplo #####
pop_haplos <- quickHaplo(nInd = 10, # number of individuals
                         nChr = 3, # number of chromosomes
                         segSites = 1000, # number of segregation sites
                         genLen = 1, # genetic length of chromosomes
                         ploidy = 2L, # ploidy level of organism
                         inbred = TRUE) # are founders inbred?

# 2. SET SIMULATION PARAMETERS #####

SP = SimParam$new(pop_haplos) # create a variable for new simulation parameters
str(SP)

#__ 2.1 Additive traits #####
SP$addTraitA(nQtlPerChr = 100, # number of QTL for trait
             mean = 500, # mean genetic value of the trait
             var = 100, # variance of trait
             corA = NULL, # matrix of correlations between additive effects
             gamma = FALSE, # to use a gamma distribution in place of a normal distribution
             shape = 1, # shape parameter for gamma distribution only
             force = FALSE) # keep false until this is understood!

# set environmental variance to create phenotypes #
SP$setVarE(h2 = 0.3, # vector of narrow sense heritability
           H2 = NULL, # vector of broad sense heritability
           varE = NULL) # vector of error variances

#__ 2.2 SNP chips #####
SP$addSnpChip(nSnpPerChr = 32, # number of SNPs per chromosome
              minSnpFreq = NULL, # minimum allowable frequency for SNP loci
              refPop = NULL) # reference population for calculating SNP frequency

# 3. DEVELOP POPULATIONS #####

gen <- 10 # specify the number of generations for which selection will be carried out

#__ 3.1 Develop population A #####
pop_A <- newPop(rawPop = pop_haplos, # object of MapPop-class
                mother = NULL, # optional id for mothers
                father = NULL, # optional id for fathers
                origM = NULL, # optional alternative id for mothers
                origF = NULL, # optional alternative id for fathers
                isDH = FALSE, # indicate if double haploids and/or inbred
                simParam = SP) # SimParam objects

POP_A <- vector(length = gen, mode = "list") # create a vector of lists to store the population for each generation
POP_A[[1]] <- pop_A # add the founder population to the vector of lists

for (i in 2:gen) {
  POP_A[[i]] <- randCross(pop = POP_A[[i - 1]],
                          nCrosses = 10,
                          nProgeny = 1,
                          balance = TRUE,
                          parents = NULL,
                          simParam = SP)}

#__ 3.2 Develop population B #####
pop_B <- newPop(rawPop = pop_haplos, # object of MapPop-class
                mother = NULL, # optional id for mothers
                father = NULL, # optional id for fathers
                origM = NULL, # optional alternative id for mothers
                origF = NULL, # optional alternative id for fathers
                isDH = FALSE, # indicate if double haploids and/or inbred
                simParam = SP) # SimParam objects

POP_B <- vector(length = gen, mode = "list") # create a vector of lists to store the population for each generation
POP_B[[1]] <- pop_B # add the founder population to the vector of lists

for (i in 2:gen) {
  POP_B[[i]] <- randCross(pop = POP_B[[i - 1]],
                          nCrosses = 10,
                          nProgeny = 1,
                          balance = TRUE,
                          parents = NULL,
                          simParam = SP)}

# 4. DIFFERENTIATE POPULATIONS #####

c_gen <- 1 # specify the generation at which PCA will be performed

snp_A <- pullSnpGeno(pop = POP_A[[c_gen]], # Pop-class object
                     snpChip = 1, # which chip to use
                     chr = NULL, # chromosome number, NULL = all
                     simParam = SP) # SimParam object

snp_B <- pullSnpGeno(pop = POP_B[[c_gen]], # Pop-class object
                     snpChip = 1, # which chip to use
                     chr = NULL, # chromosome number, NULL = all
                     simParam = SP) # SimParam object

snp_pops <- rbind(snp_A, snp_B)
snp_pops[,1:10]

pca <- prcomp(snp_pops)


pca.data <- data.frame(fish_ID = rownames(pca$x), # vector of fish IDs
                       hatchery = c(rep("A", gen),(rep("B", gen))),
                       X = pca$x[,1], # vector of PC1 in order of fish IDs
                       Y = pca$x[,2]) # vector of PC2 in order of fish IDs
 
# Assemble a data frame of PCA effects
p_comp <- 1:dim(pca$x)[2] # create a vector pf principle component IDs
pca_sd <- pca$sdev # call the standard deviation from the PCA results
pca_var <- pca_sd^2 # calculate the variance of the principle components
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1) # calculate the percent variance of the total variance for each principle component

df_p_comp <- data.frame(p_comp,
                        pca_sd,
                        pca_var,
                        pca_var_per)
head(df_p_comp)

# Create scree plot
df_p_comp %>%
  ggplot(aes(x = p_comp, y = pca_var_per)) +
  geom_point(stat = 'identity', color = "blue", size = 2) +
  geom_line(linetype = "dotted", color = "blue", size = 1) +
  theme_bw() +
  xlab("principle component") +
  ylab("effect size (percent)")


pca.data %>%
  ggplot(aes(x = X, y = Y, label = fish_ID)) +
  geom_text(size = 3) +
  geom_point(aes(color = hatchery), size = 8, alpha = 0.5) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
  scale_color_manual(values = c("A" = "blue", "B" = "grey")) +
  theme_bw() +
  ggtitle(paste("generation", c_gen, sep = " "))

# barplot(pca.var.per,
#         main="Scree Plot",
#         xlab="Principal Component",
#         ylab="Percent Variation")
# pca.data
# pca.var.per


# 5. IDENTIFY SIGNIFICANT SNPS #####

p_comp <- 1 # select the principle component to examine in more detail

effect <- TRUE # for examining POSITIVE effects

snp_id <- names(sort(pca$rotation[,p_comp],
                   decreasing = effect))
load_score <- as.numeric(sort(pca$rotation[,p_comp],
                         decreasing = effect))

df_load <- data.frame(snp_id, load_score)
df_load %>%
  mutate(ls_norm = load_score / load_score[1]) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  filter(ls_norm > 0)

effect <- FALSE # for examining NEGATIVE effects

snp_id <- names(sort(pca$rotation[,p_comp],
                     decreasing = effect))
load_score <- as.numeric(sort(pca$rotation[,p_comp],
                              decreasing = effect))

df_load <- data.frame(snp_id, load_score)
df_load %>%
  mutate(ls_norm = load_score / load_score[1]) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  filter(ls_norm > 0)

# 6. VERIFY SIGNIFICANCE OF SNPS #####

df_verify <- data.frame(hatchery = c(rep("A", 10), rep("B", 10)),
                        SNP_5 = snp_pops[,9],
                        SNP_6 = snp_pops[,17])
a <- df_verify %>%
  mutate(gType = as.character(paste(SNP_5, SNP_6, sep = "-")))
a

a$gType[1:10] %in% a$gType[11:20]
