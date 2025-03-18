library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(kableExtra)

#EEG electrode atlas cortical areas list
electrodes <- c("vlpfc", "dlpfc", "spfc", "motor")

# Delays ----
#Depth-specific (superficial/deep) R1 and fooof measures for final study sample
myelin.mse.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/R1_mse_allDelays.RDS") 
myelin.mse.7T <- lapply(myelin.mse.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_mse_electrodeatlas <- do.call(rbind, myelin.mse.7T)
SGIGR1_mse_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_mse_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 


R1.mse.maineffect.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of mse measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssocationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "sex.x"] <- "sex"

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'superficial'), mse.measure = "Var1", output.df.name = "R1Var1_superficial_maineffect.RDS")
R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'deep'), mse.measure = "Var1", output.df.name = "R1Var1_deep_maineffect.RDS")



R1.mse.interaction.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ mse*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssocationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"

R1.mse.interaction.gams(input.depth.df = SGIGR1_mse_electrodeatlas, mse.measure = "Var1", output.df.name = "R1Var1_depth_interaction.RDS")


# Plot the R1 mse effects
R1Var1_superficial_maineffect <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssocationsWithMyelin/R1Var1_superficial_maineffect.RDS")
R1Var1_superficial_maineffect$gam.statistics.df


R1Var1_deep_maineffect <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssocationsWithMyelin/R1Var1_deep_maineffect.RDS")
R1Var1_deep_maineffect$gam.statistics.df

R1mse.plotdata <- rbind(R1Var1_superficial_maineffect$gam.statistics.df %>% mutate(depth = 'superficial'), R1Var1_deep_maineffect$gam.statistics.df %>% mutate(depth = 'deep')) #long df with stats per superficial and deep depths
R1mse.plotdata$orig_parcelname <- gsub("_R1", "", R1mse.plotdata$orig_parcelname)
R1mse.plotdata <- R1mse.plotdata %>% dplyr::filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1mse.plotdata$orig_parcelname <- factor(R1mse.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1mse.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = depth)) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#90b8a6")) +
  labs(x="\n", y=sprintf("R1 - Var1 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - mse interactions by cortical depth ----

#### Interaction effects
R1Var1_depth_interaction <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssocationsWithMyelin/R1Var1_depth_interaction.RDS")

R1Var1_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 

vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_Var1, k = 3, fx = F) + s(vlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_Var1", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("VLPFC")


dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_Var1, k = 3, fx = F) + s(dlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_Var1", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("DLPFC")


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_Var1, k = 3, fx = F) + s(spfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_Var1", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("SPFC")


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_Var1, k = 3, fx = F) + s(motor_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_Var1", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("Motor")



# Rest Eyes Open----
## Depth-specific (superficial/deep) R1 and mse measures ----
myelin.mse.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/R1_mse_resteyesOpen.RDS") 
myelin.mse.7T <- lapply(myelin.mse.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_mse_electrodeatlas <- do.call(rbind, myelin.mse.7T)
SGIGR1_mse_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_mse_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.mse.maineffect.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of mse measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "sex.x"] <- "sex"

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'superficial'), mse.measure = "Var1", output.df.name = "R1Var1_superficial_eyesOpen_maineffect.RDS")

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'deep'), mse.measure = "Var1", output.df.name = "R1Var1_deep_eyesOpen_maineffect.RDS")


R1.mse.interaction.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ mse*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"

R1.mse.interaction.gams(input.depth.df = SGIGR1_mse_electrodeatlas, mse.measure = "Var1", output.df.name = "R1Var1_eyesOpen_depth_interaction.RDS")

## Plot the raw data ----
ggplot(SGIGR1_mse_electrodeatlas, aes (x = vlpfc_Var1, y = vlpfc_R1, group = depth, color = depth)) +
  geom_point(size = 1) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#90b8a6")) +
  labs(x="\n", y=sprintf("R1 - Var1 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

## Plot the R1 mse effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1Var1_superficial_eyesOpen_maineffect$gam.statistics.df
R1Var1_superficial_eyesOpen_maineffect$gam.statistics.df$depth <- 'superficial'
R1Var1_superficial_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'


R1Var1_deep_eyesOpen_maineffect$gam.statistics.df
R1Var1_deep_eyesOpen_maineffect$gam.statistics.df$depth <- 'deep'
R1Var1_deep_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'


R1mse.plotdata <- rbind(R1Var1_superficial_eyesOpen_maineffect$gam.statistics.df , R1Var1_deep_eyesOpen_maineffect$gam.statistics.df)


R1mse.plotdata$orig_parcelname <- gsub("_R1", "", R1mse.plotdata$orig_parcelname)
R1mse.plotdata <- R1mse.plotdata %>% dplyr::filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1mse.plotdata$orig_parcelname <- factor(R1mse.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1mse.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth, condition))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - Var1 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - mse interactions by cortical depth ----

#### Interaction effects

R1Var1_eyesOpen_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_Var1, k = 3, fx = F) + s(vlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_Var1[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    text = element_text(size=15),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("VLPFC")



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_Var1, k = 3, fx = F) + s(dlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_Var1 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=15),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("DLPFC")


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_Var1, k = 3, fx = F) + s(spfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("SPFC")


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_Var1, k = 3, fx = F) + s(motor_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("Motor")






# Rest Eyes Closed----
## Depth-specific (superficial/deep) R1 and mse measures ----
myelin.mse.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/R1_mse_resteyesClosed.RDS") 
myelin.mse.7T <- lapply(myelin.mse.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})


SGIGR1_mse_electrodeatlas <- do.call(rbind, myelin.mse.7T)
SGIGR1_mse_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_mse_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.mse.maineffect.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of mse measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "sex.x"] <- "sex"

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'superficial'), mse.measure = "Var1", output.df.name = "R1Var1_superficial_eyesClosed_maineffect.RDS")

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'deep'), mse.measure = "Var1", output.df.name = "R1Var1_deep_eyesClosed_maineffect.RDS")


R1.mse.interaction.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ mse*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"

R1.mse.interaction.gams(input.depth.df = SGIGR1_mse_electrodeatlas, mse.measure = "Var1", output.df.name = "R1Var1_eyesClosed_depth_interaction.RDS")





## Plot the R1 mse effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1Var1_superficial_eyesClosed_maineffect$gam.statistics.df
R1Var1_superficial_eyesClosed_maineffect$gam.statistics.df$depth <- 'superficial'
R1Var1_superficial_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'


R1Var1_deep_eyesClosed_maineffect$gam.statistics.df
R1Var1_deep_eyesClosed_maineffect$gam.statistics.df$depth <- 'deep'
R1Var1_deep_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'


R1mse.plotdata <- rbind(R1Var1_superficial_eyesClosed_maineffect$gam.statistics.df , R1Var1_deep_eyesClosed_maineffect$gam.statistics.df)

R1mse.plotdata$orig_parcelname <- gsub("_R1", "", R1mse.plotdata$orig_parcelname)
R1mse.plotdata <- R1mse.plotdata %>% dplyr::filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1mse.plotdata$orig_parcelname <- factor(R1mse.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1mse.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth, condition))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - Var1 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - mse interactions by cortical depth ----

#### Interaction effects

R1Var1_eyesClosed_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_Var1, k = 3, fx = F) + s(vlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_Var1[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    text = element_text(size=15),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("VLPFC")



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_Var1, k = 3, fx = F) + s(dlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_Var1 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=15),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("DLPFC")


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_Var1, k = 3, fx = F) + s(spfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("SPFC")


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_Var1, k = 3, fx = F) + s(motor_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("Motor")



# Fix ----
## Depth-specific (superficial/deep) R1 and mse measures ----
myelin.mse.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/R1_mse_fix.RDS") 
myelin.mse.7T <- lapply(myelin.mse.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_mse_electrodeatlas <- do.call(rbind, myelin.mse.7T)
SGIGR1_mse_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_mse_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.mse.maineffect.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of mse measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- mse.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "sex.x"] <- "sex"

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'superficial'), mse.measure = "Var1", output.df.name = "R1Var1_superficial_fix_maineffect.RDS")

R1.mse.maineffect.gams(input.depth.df = SGIGR1_mse_electrodeatlas %>% dplyr::filter(depth == 'deep'), mse.measure = "Var1", output.df.name = "R1Var1_deep_fix_maineffect.RDS")


R1.mse.interaction.gams <- function(input.depth.df, mse.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ mse*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, mse.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_mse_electrodeatlas)[colnames(SGIGR1_mse_electrodeatlas) == "age.x"] <- "age"

R1.mse.interaction.gams(input.depth.df = SGIGR1_mse_electrodeatlas, mse.measure = "Var1", output.df.name = "R1Var1_fix_depth_interaction.RDS")



## Plot the R1 mse effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1Var1_superficial_fix_maineffect$gam.statistics.df
R1Var1_superficial_fix_maineffect$gam.statistics.df$depth <- 'superficial'


R1Var1_deep_fix_maineffect$gam.statistics.df
R1Var1_deep_fix_maineffect$gam.statistics.df$depth <- 'deep'


R1mse.plotdata <- rbind(R1Var1_superficial_fix_maineffect$gam.statistics.df , R1Var1_deep_fix_maineffect$gam.statistics.df)


R1mse.plotdata$orig_parcelname <- gsub("_R1", "", R1mse.plotdata$orig_parcelname)
R1mse.plotdata <- R1mse.plotdata %>% dplyr::filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1mse.plotdata$orig_parcelname <- factor(R1mse.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1mse.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - Var1 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - mse interactions by cortical depth ----

#### Interaction effects

R1Var1_fix_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_Var1, k = 3, fx = F) + s(vlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_Var1[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    text = element_text(size=15),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("VLPFC")



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_Var1, k = 3, fx = F) + s(dlpfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_Var1 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=15),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("DLPFC")


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_Var1, k = 3, fx = F) + s(spfc_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("SPFC")


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_Var1, k = 3, fx = F) + s(motor_Var1, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_mse_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_Var1[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nVar1", y=sprintf("R1 (adjusted)\n")) +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 15, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 15, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    plot.title = element_text(hjust = 0.5))+ ggtitle("Motor")

