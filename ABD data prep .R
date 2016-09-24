# ------------------------------------------------------------------------------
# Program: Data Prep and Non-twin Analyses for ABD sample
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
  

# Graphical display ----------------------------------------------------------------------------------------------------------------------------------------------------------  
  
# participation indicators
  dat$wave1 <- 1 # W1
  dat$wave2 <- ifelse(is.na(dat$W2AGE), 0, 1) # W2
  dat$wave3 <- ifelse(is.na(dat$w3age), 0, 1) # W3
  dat$wave4 <- ifelse(is.na(dat$w4age), 0, 1) # W4
  
# Alcohol initiation: alce8 - alce17
  init_vars <- c(rep(paste0("alce",8:17)))
# Intoxication frequency
  intox_vars <- c(rep(paste0("toxe",8:17)))
# Use in past 90 days: : alcc8 - alcc17
  use_vars <- c(rep(paste0("alcc",8:17)))
# Alcohol problems:
  #Checking totals
  t1 <- table(dat$aralcsx1, dat$W1AGE)
  t2 <- table(dat$aralcsx2, dat$W2AGE)
  t2 <- cbind(matrix(c(rep(0,2)), nrow = 2, ncol = 1), t2)
  t3 <- table(dat$aralcsx3, dat$w3age)
  t3 <- cbind(matrix(c(rep(0,8)), nrow = 2, ncol = 4), t3, c(0,0))
  t4 <- table(dat$aralcsx4, dat$w4age)
  t4 <- cbind(matrix(c(rep(0,12)), nrow = 2, ncol = 6), t4, c(0,0))
  (checks <- rbind(t1, t2, t3, t4))
  probs <- checks[c(2,4,6,8),]
  noprobs <- checks[c(1,3,5,7),]
  rbind(colSums(probs), colSums(noprobs))
  
  prob_vars <- c(rep(paste0("alcprob",8:17)))
    for (i in 8:17){
        dat[, paste0("alcprob",i)] <- NA
      }
  age_vars  <-c("W1AGE", "W2AGE", "w3age", "w4age")
    for (i in 1:4){
      dat[, age_vars[i]][is.na(dat[, age_vars[i]])] <- -9
    }
  wave_vars <-c(rep(paste0("wave",1:4)))
  adsx_vars <-c(rep(paste0("aralcsx",1:4)))

  # Loops for restructuring on age rather than wave (filling values into alcprob8 - alcprob17)  
    for (j in 1:4){  
      for (i in 8:17) { 
        dat[, paste0("alcprob",i)] <- ifelse(dat[, age_vars[j]] == i & !is.na(dat[, adsx_vars[j]]), 
                                             dat[,adsx_vars[j]], 
                                             dat[, paste0("alcprob",i)])  
      }
    }
  
  
  lapply(dat[, prob_vars], table)


  
# Combined 
  alc_vars <- c(init_vars, use_vars, intox_vars, prob_vars)
  
# Calculating proportion of endorsement across each age for each alcohol phenotype
mean <- NULL    
  for (i in 1:length(alc_vars)){
      x <- round(mean(dat[, alc_vars[i]], na.rm = T), 3)
      mean <- rbind(mean, x)
    }
age <- as.numeric(rep(c(8:17), 4)) 
phenotype <- c(rep("Initiation", 10), rep("Use (past 90 days", 10), rep("Intoxication", 10), rep("Problems", 10))
alc_means <- as.data.frame(cbind(phenotype, age, mean), row.names = F, stringsAsFactors = F)

ggplot(alc_means, aes(as.numeric(age), as.numeric(V3)))  + 
  geom_point(aes(colour = phenotype, shape = phenotype)) +
  geom_line(aes(group = phenotype, colour = phenotype)) + 
    expand_limits(y=c(0,.75), x=c(8:17)) + 
    xlab('Age') + 
    ylab('% Endorsing') + 
    ggtitle('Alcohol Related Phenotypes from Ages 8 to 17')
    

# Clearing up the work space
rm(alc_means, checks, mean, noprobs, probs, t1, t2, t3, t4, age, age_vars, i, j, phenotype, wave_vars, adsx_vars, x)  
  


# ML growth curve models ----------------------------------------------------------------------------------------------------------------------------

# Reshaping data
gc_dat <- subset(dat, select = c("FAMNO", 'twid', 'poor', alc_vars))
gc_dat$ID <- paste0(gc_dat$FAMNO, gc_dat$twid)
gc_dat$twid <- NULL

# Subsetting, reshaping, and merging across phenotype using reshape2 (HW)
gc_init <- subset(gc_dat, select = c("FAMNO", 'ID', 'poor', init_vars))
  gc_init <- melt(gc_init, id.vars = c("FAMNO", "ID", 'poor'), measured = init_vars, variable.name = "age", value.name = "initiation")
  gc_init$age <- as.numeric(substr(gc_init$age, 5, 6)) # substringing age to numeric value  
  
gc_use <- subset(gc_dat, select = c("FAMNO", 'ID', 'poor', use_vars))
  gc_use <- melt(gc_use, id.vars = c("FAMNO", "ID", 'poor'), measured = use_vars, variable.name = "age", value.name = "use")
  gc_use$age <- as.numeric(substr(gc_use$age, 5, 6)) # substringing age to numeric value  
  
gc_intx <- subset(gc_dat, select = c("FAMNO", 'ID', 'poor', intox_vars))
  gc_intx <- melt(gc_intx, id.vars = c("FAMNO", "ID", 'poor'), measured = intox_vars, variable.name = "age", value.name = "intox")
  gc_intx$age <- as.numeric(substr(gc_intx$age, 5, 6)) # substringing age to numeric value  
  
gc_prob <- subset(gc_dat, select = c("FAMNO", 'ID', 'poor', prob_vars))
  gc_prob <- melt(gc_prob, id.vars = c("FAMNO", "ID", 'poor'), measured = prob_vars, variable.name = "age", value.name = "problems")
  gc_prob$age <- as.numeric(substr(gc_prob$age, 8, 9)) # substringing age to numeric value  


# Merging across phenotype for a single long dataset.  
gc_dat_long <- Reduce(function(x, y) merge(x, y, by = c("FAMNO", 'ID', 'poor', 'age'), all=TRUE), list(gc_init, gc_use, gc_intx, gc_prob))
rm(gc_init, gc_use, gc_intx, gc_prob, gc_dat) # cleaning up the global env.


#Unconditional model
glmm1 <- glmer(intox ~ 1 + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm1)

glmm2 <- glmer(intox ~ as.numeric(age) + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm2)

glmm3 <- glmer(intox ~ as.numeric(age) + (as.numeric(age) | ID), family = "binomial", data = gc_dat_long)
summary(glmm3)
