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
        dat[, paste0("alcprob",i)] <- -9
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
  alc_vars <- c(init_vars, use_vars, intox_vars)
  
# Calculating proportion of endorsement across each age for each alcohol phenotype
mean <- NULL    
  for (i in 1:length(alc_vars)){
      x <- round(mean(dat[, alc_vars[i]], na.rm = T), 3)
      mean <- rbind(mean, x)
    }
age <- as.numeric(rep(c(8:17), 3)) 
pheno <- c(rep("Init", 10), rep("Use", 10), rep("Intox", 10))
alc_means <- as.data.frame(cbind(pheno, age, mean), row.names = F, stringsAsFactors = F)

ggplot(alc_means, aes(as.numeric(age), as.numeric(V3)))  + 
  geom_point(aes(colour = pheno, shape = pheno)) +
  geom_line(aes(group = pheno, colour = pheno)) + 
    expand_limits(y=c(0,1), x=c(8:17)) + 
    xlab('Age') + 
    ylab('% Endorsing') 
    
  
  

# ML growth curve models ----------------------------------------------------------------------------------------------------------------------------

# Reshaping data
gc_dat <- subset(dat, select = c("FAMNO", 'twid', 'poor', intox_vars[1:10]))
gc_dat$ID <- paste0(gc_dat$FAMNO, gc_dat$twid)
gc_dat$twid <- NULL

# Using reshape
# gc_dat_long <- reshape(gc_dat, idvar = c("FAMNO", 'ID'), varying = intox_vars[1:10], v.names = "intox", timevar = "age", times = as.character(c(8:17)), direction = "long")

#Using reshape2 (HW)
gc_dat_long <- melt(gc_dat, id.vars = c("FAMNO", "ID", 'poor'), measured = intox_vars[1:10], variable.name = "age", value.name = "intox")
gc_dat_long$age <- substr(gc_dat_long$age, 5, 6)

#Unconditional model
glmm1 <- glmer(intox ~ 1 + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm1)

glmm2 <- glmer(intox ~ as.numeric(age) + (1 | ID), family = "binomial", data = gc_dat_long, nAGQ = 10)
summary(glmm2)

glmm3 <- glmer(intox ~ as.numeric(age) + (as.numeric(age) | ID), family = "binomial", data = gc_dat_long)
summary(glmm3)
