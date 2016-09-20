# ------------------------------------------------------------------------------
# Program: Data Prep and Non-twin Analyses for ABD sample
#  Author: Peter Barr
#    Date: 9 18 2016
#
# -------|---------|---------|---------|---------|---------|---------|---------|

rm(list= ls())
setwd("~/Desktop/R/ABD/")

# R packages
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


  
# Alcohol initiation: alce8 - alce17, alcey, alcez
init_vars <- grep(x = names(dat), pattern = "^alce[0-9A-z]", value = TRUE)
  init_vars <- init_vars[7:18]

use_vars <- grep(x = names(dat), pattern = "^alcc[0-9]", value = TRUE)

days_dr <- grep(x = names(dat), pattern = "^alcnd[0-9]", value = TRUE)

intox_vars <- grep(x = names(dat), pattern = "^toxe[0-9A-z]", value = TRUE)
intox_vars <- intox_vars[-c(1,2,15:19)]


# Alcohol inititation
plots <- list()  # new empty list
for (i in 1:10) {
  p <- qplot(dat[,init_vars[i]], main = init_vars[i])
  plots[[i]] <- p  # add each plot into plot list
}
multiplot(plotlist = plots, cols = 2)






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
