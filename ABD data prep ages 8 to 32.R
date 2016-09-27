# ------------------------------------------------------------------------------
# Program: Data Prep and Non-twin Analyses for ABD sample: Analyses for ages 8 -32
#  Author: Peter Barr
#    Date: 9 18 2016
#
# -------|---------|---------|---------|---------|---------|---------|---------|

rm(list= ls())
setwd("~/Desktop/R/ABD/")

require(psych)
require(polycor)
require(ggplot2)
require(lme4)
require(reshape2)
require(sas7bdat)
source("~/Desktop/R/Functions/functions.R")

# Read in ABD data ----------------------------------------------------------------------------------------------------------------------------------------------------------

# ABD data
dat <- read.csv("data/ABDphens9.csv", header = TRUE)
# Existing neighborhood data
census <- read.csv("data/abd_census.csv", header = TRUE)
# Merge data
dat <- merge(dat, census, by = "FAMNO", all.x = TRUE)
variable.names(dat)
rm(census)

# dropping cases with missing age at Wave I (lose 28 cases)
dat <- subset(dat, !is.na(W1AGE))

# participation indicators
dat$wave1 <- 1 # W1
dat$wave2 <- ifelse(is.na(dat$W2AGE), 0, 1) # W2
dat$wave3 <- ifelse(is.na(dat$w3age), 0, 1) # W3
dat$wave4 <- ifelse(is.na(dat$w4age), 0, 1) # W4
dat$yafu[is.na(dat$yafu)] <- 0
dat$tsa[is.na(dat$tsa)] <- 0

dat$participation <- dat$wave1 + dat$wave2 + dat$wave3 + dat$wave4 + dat$yafu + dat$tsa 
table(dat$participation)


# Graphical display ----------------------------------------------------------------------------------------------------------------------------------------------------------  

# Initiation and problems are the only variables with consistent measures across all six time points

alcinit_vars <- c(rep(paste0("aralce",1:4)), 'alcevery', 'alceverz') # initiation
adsx_vars <-c(rep(paste0("aralcnr",1:4)), 'alcnimy', 'alcnimz') # problems
for (j in 1:6){
  dat[, adsx_vars[j]] <- ifelse(dat[, alcinit_vars[j]] == 0 & is.na(dat[, adsx_vars[j]]), 0, dat[, adsx_vars[j]])
}
age_vars  <-c("W1AGE", "W2AGE", "w3age", "w4age", "wyage", 'wzage') # age at each wave


# REcoding missing age as -9
  for (i in 1:6){
    dat[, age_vars[i]][is.na(dat[, age_vars[i]])] <- -9
  }

# Alcohol initiation: alcinit8 - alcinit32
  init_vars <- c(rep(paste0("alcinit",8:32)))
  for (i in 8:32){
    dat[, paste0("alcinit",i)] <- NA
  }
  
  # Loops for restructuring on age rather than wave (filling values into alcprob8 - alcprob32)  
  for (j in 1:6){  
    for (i in 8:32) { 
      dat[, paste0("alcinit",i)] <- ifelse(dat[, age_vars[j]] == i & !is.na(dat[, alcinit_vars[j]]), 
                                           dat[,alcinit_vars[j]], 
                                           dat[, paste0("alcinit",i)])
    }
  }
  
  lapply(dat[, init_vars], table)
  
  
# Alcohol problems:
  prob_vars <- c(rep(paste0("alcprob",8:32)))
  for (i in 8:32){
    dat[, paste0("alcprob",i)] <- NA
  }
  
# Loops for restructuring on age rather than wave (filling values into alcprob8 - alcprob32)  
  for (j in 1:6){  
    for (i in 8:32) { 
      dat[, paste0("alcprob",i)] <- ifelse(dat[, age_vars[j]] == i & !is.na(dat[, adsx_vars[j]]), 
                                           dat[,adsx_vars[j]], 
                                           dat[, paste0("alcprob",i)])
      dat[, paste0("alcprob",i)] <- ifelse(dat[, paste0("alcprob",i)] > 7, 7, dat[, paste0("alcprob",i)])  
    }
  }
  
  lapply(dat[, prob_vars], table)
  

  
# Calculating proportion of endorsement across each age for each alcohol phenotype
alc_vars <- c(init_vars, prob_vars)
mean <- NULL    
for (i in 1:length(alc_vars)){
  x <- round(mean(dat[, alc_vars[i]], na.rm = T), 3)
  mean <- rbind(mean, x)
}
age <- c(rep(8:32, 2))
phenotype <- c(rep("Initiation", 25), rep("Problems", 25))
alc_means <- as.data.frame(cbind(phenotype, age, mean), row.names = F, stringsAsFactors = F)

ggplot(alc_means, aes(x = as.numeric(age), y =as.numeric(V3), colour = phenotype))  + 
  geom_point() +
  geom_line() +
  expand_limits(y=c(0,1.5)) + 
  xlab('Age') + 
  ylab('Mean Problems') 




# Clearing up the work space
rm(alc_means, mean, age, age_vars, i, j, phenotype, alcinit_vars, adsx_vars, x)  




# ML growth curve models ----------------------------------------------------------------------------------------------------------------------------

# Covariates
cov_vars <- c("pmon", "ppcon", "pcrel", "pmddxi", "pgadxi", "palcxi", "pintxi", "pextxi", "pallxi","matage",
              "magelt18", "popoc", "momoc", "poplowoc", "momlowoc", "poped", "momed", "poplowed", "momlowed", 
              "income", "ksib", "sesyear", "FATQ11", "singleparent", "famsize", "incomelevel", "povertylevel",
              "poor", "t1bwt", "t2bwt", "t1lowbwt", "t2lowbwt", "pregdrink", "pregsmoke", "URBAN_1", "RURAL_1",
              "HIGH1", "COLLEGE1", "POV", "T_WORK", "MINORITY", "HOUSE", "NWORK_M", "UNEMPL_M", "NWORK_W",
              "UNEMPL_W", "PERCAP", "WPERCAP", "INC_H", "INC_F")

# Reshaping data
gc_dat <- subset(dat, select = c("FAMNO", 'twid', cov_vars, alc_vars))
gc_dat$ID <- paste0(gc_dat$FAMNO, gc_dat$twid)
gc_dat$twid <- NULL

# Subsetting, reshaping, and merging across phenotype using reshape2 (HW)
gc_init <- subset(gc_dat, select = c("FAMNO", 'ID', cov_vars, init_vars))
gc_init <- melt(gc_init, id.vars = c("FAMNO", "ID", cov_vars), measured = init_vars, variable.name = "age", value.name = "initiation")
gc_init$age <- as.numeric(substr(gc_init$age, 8, 9)) # substringing age to numeric value  

gc_prob <- subset(gc_dat, select = c("FAMNO", 'ID', cov_vars, prob_vars))
gc_prob <- melt(gc_prob, id.vars = c("FAMNO", "ID", cov_vars), measured = prob_vars, variable.name = "age", value.name = "problems")
gc_prob$age <- as.numeric(substr(gc_prob$age, 8, 9)) # substringing age to numeric value  


# Merging across phenotype for a single long dataset.  
gc_dat_long <- Reduce(function(x, y) merge(x, y, by = c("FAMNO", 'ID', cov_vars, 'age'), all=TRUE), list(gc_init, gc_prob))
rm(gc_init, gc_prob, gc_dat) # cleaning up the global env.


# Saving data for import into STATA
write.csv(gc_dat_long, "~/Desktop/ABD alcohol phenotypes long.csv", row.names = F, na = ".")


#Growth models for alcohol use inititation
glmm1 <- glmer(initiation ~ 1 + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm1)

glmm2 <- glmer(initiation ~ age + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm2)

glmm3 <- glmer(initiation ~ age + (age | ID), family = binomial(link = "logit"), data = gc_dat_long)
summary(glmm3)


