# load libraries

library(arm)
library(readxl)
library(Matrix)
library(lme4)
library(car)
library(ggplot2)
library(dplyr) # for pipe
library(tidyr)
library(lmerTest)
library(sjPlot)
library(MASS) # for rlm
library(caret)
library(robustbase) 
library(svglite)
library(raincloudplots)
library(tidyverse) 
library(performance) # to check model assumptions
library(see)
library(emmeans)
library(effects)
library(bayestestR)
library(brms)
library(BAS)
library(insight)


## load in data
tracking1 <- read.csv("tracking_all.csv")
tracking1 <- tracking1[, c("subID", "read.enjoyment", "Age", "read.time")]
tracking1 <- tracking1[-c(19, 26, 31), ] # remove outliers identified in previous steps
median <- read.csv("behaviouralclean_RTmedian.csv") # load RT data and GMSI

# load MI data
MIdata <- read.csv("behclus_RT_newlags.txt")
MIdata <- MIdata[ ,-1]
median <- median[ ,-1]

tracking1 <-  cbind(tracking1, median, MIdata)

# change names for pivoting
names(tracking1)[names(tracking1) == "MIspeechneg"] <- "MI_speech"
names(tracking1)[names(tracking1) == "MImusicneg"] <- "MI_music"

# rearrange data to long format
tracking_long <- tracking1 %>%
  pivot_longer(
    cols = starts_with("MI_"),
    names_to = c(".value", "stimulus"),
    names_sep = "_"
  )

# normalise
vars <- c(
  "MI",
  "Age", 
  "GMSI", 
  "RTmedianall",
  "read.enjoyment",
  "read.time")

tracking_long <- tracking_long %>%
  mutate(across(all_of(vars), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

tracking_long$stimulus <- relevel(factor(tracking_long$stimulus), ref = "speech")

# model we are using -> model 1
model1 <- rlm(RTmedianall_z ~ MI_z*stimulus*GMSI_z + Age_z + read.time_z+read.enjoyment_z, data = tracking_long)

# model without interaction
model2 <- rlm(RTmedianall_z ~ MI_z*GMSI_z + stimulus*GMSI_z + Age_z + read.enjoyment_z + read.time_z, data = tracking_long)

# bayesian factor
bayesfactor_models(model1,denominator = model2, verbose = TRUE)


