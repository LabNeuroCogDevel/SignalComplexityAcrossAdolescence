

# Load Libraries 
library(LNCDR)
library(dplyr)
library(ggplot2)
library(mgcv)
library(lme4)
library(lmerTest)
library(tidyr)
library(gratia)
library(tidyverse)
library(ggeffects)
library(marginaleffects)

behavioralInteraction <- function(df, location, brainMeasure, behavioralMeasure){
  df <- df %>% dplyr::filter(region == location)
  df$epoch <- factor(df$epoch, ordered=T)

  gamModelformula <- as.formula(sprintf("%s ~ s(age, k = 3, fx = F) + %s*epoch + epoch + sex + s(lunaid, bs='re')", behavioralMeasure, brainMeasure))

  model <- gam(gamModelformula, method = c("REML"), data = df) 
  
  summary_model <- summary(model)
  df_results <- as.data.frame(summary_model$p.table)  # Extract parametric coefficients
  df_smooth <- as.data.frame(summary_model$s.table)   # Extract smooth terms

  resultList <- list(model, df_results, df_smooth)
  names(resultList) <- list("model", "fixedResults", "smoothedResults")
  
  return(resultList)
} 

behavioralAgeInteraction <- function(df, location, brainMeasure, behavioralMeasure){
  df <- df %>% dplyr::filter(region == location)
  df$epoch <- factor(df$epoch, ordered=T)
  
  gamModelformula <- as.formula(sprintf("%s ~ s(age, k = 3, fx = F) + s(age, k = 3, fx = F, by = %s) + epoch + sex + s(lunaid, bs='re')", behavioralMeasure, brainMeasure))
  
  model <- gam(gamModelformula, method = c("REML"), data = df) 
  
  summary_model <- summary(model)
  df_results <- as.data.frame(summary_model$p.table)  # Extract parametric coefficients
  df_smooth <- as.data.frame(summary_model$s.table)   # Extract smooth terms
  
  resultList <- list(model, df_results, df_smooth)
  names(resultList) <- list("model", "fixedResults", "smoothedResults")
  
  return(resultList)
} 

behavioralRegionInteraction <- function(df, location, brainMeasure, behavioralMeasure){
  df$region <- factor(df$region, ordered=T)
  
  gamModelformula <- as.formula(sprintf("%s ~ s(age, k = 3, fx = F) + %s*region + region + epoch + sex + s(lunaid, bs='re')", behavioralMeasure, brainMeasure))
  
  model <- gam(gamModelformula, method = c("REML"), data = df) 
  
  summary_model <- summary(model)
  df_results <- as.data.frame(summary_model$p.table)  # Extract parametric coefficients
  df_smooth <- as.data.frame(summary_model$s.table)   # Extract smooth terms
  
  resultList <- list(model, df_results, df_smooth)
  names(resultList) <- list("model", "fixedResults", "smoothedResults")
  
  return(resultList)
} 

behavioralMainEffect <- function(df, location, brainMeasure, behavioralMeasure){
  df <- df %>% dplyr::filter(region == location)

  gamModelformula <- as.formula(sprintf("%s ~ s(age, k = 4, fx = F) + %s + epoch + sex + s(lunaid, bs='re')", behavioralMeasure, brainMeasure))
  
  model <- gam(gamModelformula, method = c("REML"), data = df) 
  
  summary_model <- summary(model)
  df_results <- as.data.frame(summary_model$p.table)  # Extract parametric coefficients
  df_smooth <- as.data.frame(summary_model$s.table)   # Extract smooth terms
  
  resultList <- list(model, df_results, df_smooth)
  names(resultList) <- list("model", "fixedResults", "smoothedResults")
  
  return(resultList)
} 
