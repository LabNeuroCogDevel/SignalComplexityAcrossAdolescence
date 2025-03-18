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
myelin.acw.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/R1_ACW_allDelays.RDS") 
myelin.acw.7T <- lapply(myelin.acw.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_acw_electrodeatlas <- do.call(rbind, myelin.acw.7T)
SGIGR1_acw_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_acw_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 


R1.acw.maineffect.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of acw measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssocationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "sex.x"] <- "sex"

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_maineffect.RDS")
R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_maineffect.RDS")




R1.acw.interaction.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ acw*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssocationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"

R1.acw.interaction.gams(input.depth.df = SGIGR1_acw_electrodeatlas, acw.measure = "ACW_50", output.df.name = "R1ACW50_depth_interaction.RDS")
 

# Plot the R1 ACW effects
R1ACW_50_superficial_maineffect <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssocationsWithMyelin/R1ACW_50_superficial_maineffect.RDS")
R1ACW_50_superficial_maineffect$gam.statistics.df


R1ACW_50_deep_maineffect <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssocationsWithMyelin/R1ACW_50_deep_maineffect.RDS")
R1ACW_50_deep_maineffect$gam.statistics.df

R1acw.plotdata <- rbind(R1ACW_50_superficial_maineffect$gam.statistics.df %>% mutate(depth = 'superficial'), R1ACW_50_deep_maineffect$gam.statistics.df %>% mutate(depth = 'deep')) #long df with stats per superficial and deep depths
R1acw.plotdata$orig_parcelname <- gsub("_R1", "", R1acw.plotdata$orig_parcelname)
R1acw.plotdata <- R1acw.plotdata %>% filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1acw.plotdata$orig_parcelname <- factor(R1acw.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1acw.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = depth)) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#90b8a6")) +
  labs(x="\n", y=sprintf("R1 - ACW_50 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - acw interactions by cortical depth ----

#### Interaction effects
R1ACW50_depth_interaction <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssocationsWithMyelin/R1ACW50_depth_interaction.RDS")

R1ACW50_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 

vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_ACW_50, k = 3, fx = F) + s(vlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_ACW_50", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_ACW_50, k = 3, fx = F) + s(dlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_ACW_50", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_ACW_50, k = 3, fx = F) + s(spfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_ACW_50", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_ACW_50, k = 3, fx = F) + s(motor_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_ACW_50", "depth"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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
## Depth-specific (superficial/deep) R1 and acw measures ----
myelin.acw.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/R1_ACW_resteyesOpen.RDS") 
myelin.acw.7T <- lapply(myelin.acw.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_acw_electrodeatlas <- do.call(rbind, myelin.acw.7T)
SGIGR1_acw_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_acw_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.acw.maineffect.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of acw measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "sex.x"] <- "sex"

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_eyesOpen_maineffect.RDS")

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_eyesOpen_maineffect.RDS")


R1.acw.interaction.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ acw*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"

R1.acw.interaction.gams(input.depth.df = SGIGR1_acw_electrodeatlas, acw.measure = "ACW_50", output.df.name = "R1ACW50_eyesOpen_depth_interaction.RDS")

## Plot the raw data ----
ggplot(SGIGR1_acw_electrodeatlas, aes (x = vlpfc_ACW_50, y = vlpfc_R1, group = depth, color = depth)) +
  geom_point(size = 1) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#90b8a6")) +
  labs(x="\n", y=sprintf("R1 - ACW_50 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

## Plot the R1 ACW effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df
R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df$depth <- 'superficial'
R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'

R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df
R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df$depth <- 'superficial'
R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'

R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df
R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df$depth <- 'deep'
R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'

R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df
R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df$depth <- 'deep'
R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'


R1acw.plotdata <- rbind(R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df , R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df) %>%
  rbind(.,R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df) %>% 
  rbind(., R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df)


R1acw.plotdata$orig_parcelname <- gsub("_R1", "", R1acw.plotdata$orig_parcelname)
R1acw.plotdata <- R1acw.plotdata %>% filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1acw.plotdata$orig_parcelname <- factor(R1acw.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1acw.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth, condition))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - ACW_50 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - acw interactions by cortical depth ----

#### Interaction effects

R1ACW50_eyesOpen_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_ACW_50, k = 3, fx = F) + s(vlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_ACW_50[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_ACW_50, k = 3, fx = F) + s(dlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_ACW_50 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_ACW_50, k = 3, fx = F) + s(spfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_ACW_50, k = 3, fx = F) + s(motor_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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
## Depth-specific (superficial/deep) R1 and acw measures ----
myelin.acw.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/R1_ACW_resteyesClosed.RDS") 
myelin.acw.7T <- lapply(myelin.acw.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_acw_electrodeatlas <- do.call(rbind, myelin.acw.7T)
SGIGR1_acw_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_acw_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.acw.maineffect.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of acw measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "sex.x"] <- "sex"

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_eyesClosed_maineffect.RDS")

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_eyesClosed_maineffect.RDS")


R1.acw.interaction.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ acw*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"

R1.acw.interaction.gams(input.depth.df = SGIGR1_acw_electrodeatlas, acw.measure = "ACW_50", output.df.name = "R1ACW50_eyesClosed_depth_interaction.RDS")







## Plot the R1 ACW effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df
R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df$depth <- 'superficial'
R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'

R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df
R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df$depth <- 'superficial'
R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'

R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df
R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df$depth <- 'deep'
R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df$condition <- 'eyesOpen'

R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df
R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df$depth <- 'deep'
R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df$condition <- 'eyesClosed'


R1acw.plotdata <- rbind(R1ACW_50_superficial_eyesOpen_maineffect$gam.statistics.df , R1ACW_50_deep_eyesOpen_maineffect$gam.statistics.df) %>%
  rbind(.,R1ACW_50_superficial_eyesClosed_maineffect$gam.statistics.df) %>% 
  rbind(., R1ACW_50_deep_eyesClosed_maineffect$gam.statistics.df)


R1acw.plotdata$orig_parcelname <- gsub("_R1", "", R1acw.plotdata$orig_parcelname)
R1acw.plotdata <- R1acw.plotdata %>% filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1acw.plotdata$orig_parcelname <- factor(R1acw.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1acw.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth, condition))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - ACW_50 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - acw interactions by cortical depth ----

#### Interaction effects

R1ACW50_eyesClosed_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_ACW_50, k = 3, fx = F) + s(vlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_ACW_50[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_ACW_50, k = 3, fx = F) + s(dlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_ACW_50 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_ACW_50, k = 3, fx = F) + s(spfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_ACW_50, k = 3, fx = F) + s(motor_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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
## Depth-specific (superficial/deep) R1 and acw measures ----
myelin.acw.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/R1_ACW_fix.RDS") 
myelin.acw.7T <- lapply(myelin.acw.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_acw_electrodeatlas <- do.call(rbind, myelin.acw.7T)
SGIGR1_acw_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_acw_electrodeatlas)), levels = c("deep", "superficial"), ordered = T) 

R1.acw.maineffect.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.covariatesmooth.maineffect function to get statistics for a main effect of acw measure on R1 as well as fitted values/smooth estimates describing this relationship
  gam.outputs.regionlist <- lapply(electrodes, function(r){  #list of gam results, list elements are cortical areas
    gam.covariatesmooth.maineffect(input.df = input.depth.df, region = sprintf("%s_R1", r), smooth_var = 'age', smooth_var_knots = 4, 
                                   smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                   linear_covariates = "NA", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  #Extract and combine gam statistics ouputs
  gam.statistics.df <- lapply(gam.outputs.regionlist, '[[', "gam.statistics" ) #extract this df from each region's list
  gam.statistics.df <- do.call(rbind, gam.statistics.df) #merge them into one 
  
  #Extract and combine fitted values 
  gam.fittedvalues.df <- lapply(gam.outputs.regionlist, '[[', "gam.fittedvalues" ) #extract this df from each region's list
  gam.fittedvalues.df <- lapply(gam.fittedvalues.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.fittedvalues.df <- do.call(rbind, gam.fittedvalues.df) #merge them into one 
  
  #Extract and combine smooth estimates
  gam.smoothestimates.df <- lapply(gam.outputs.regionlist, '[[', "gam.smoothestimates" ) #extract this df from each region's list
  gam.smoothestimates.df <- lapply(gam.smoothestimates.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.smoothestimates.df <- do.call(rbind, gam.smoothestimates.df) #merge them into one 
  
  #Extract and combine derivatives
  gam.derivatives.df <- lapply(gam.outputs.regionlist, '[[', "gam.derivatives" ) #extract this df from each region's list
  gam.derivatives.df <- lapply(gam.derivatives.df, function(output){
    region.name <- sub("_.*", "", names(output)[1]) 
    output$region <- region.name
    names(output)[1] <- acw.measure
    return(output)
  })
  gam.derivatives.df <- do.call(rbind, gam.derivatives.df) #merge them into one 
  
  #Save depth-specific results as an RDS
  gam.outputs.statslist <- list(gam.statistics.df, gam.fittedvalues.df, gam.smoothestimates.df, gam.derivatives.df) #list of gam results, list elements are stats dfs binded across regions
  names(gam.outputs.statslist) <- list("gam.statistics.df", "gam.fittedvalues.df", "gam.smoothestimates.df", "gam.derivatives.df")
  saveRDS(gam.outputs.statslist, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
  
  #Clean up
  gc(gam.outputs.statslist)
}


colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"
colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "sex.x"] <- "sex"

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_fix_maineffect.RDS")

R1.acw.maineffect.gams(input.depth.df = SGIGR1_acw_electrodeatlas %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_fix_maineffect.RDS")


R1.acw.interaction.gams <- function(input.depth.df, acw.measure, output.df.name){
  
  #Run the gam.factorsmooth.interaction function to get R1 ~ acw*depth interaction statistics
  
  gam.outputs.regionlist <- lapply(electrodes, function(r){ 
    gam.factorsmooth.interaction(input.df = input.depth.df, region = sprintf("%s_R1", r), 
                                 smooth_var = "age", smooth_var_knots = 4, smooth_covariate = sprintf("%s_%s", r, acw.measure), smooth_covariate_knots = 3, 
                                 int_var = "depth", linear_covariates = "depth", id_var = "subject_id", random_intercepts = TRUE, set_fx = FALSE)}) 
  
  
  #Extract and combine base effect outputs
  gam.baseeffects.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.baseeffect" ) #extract this df from each region's list
  gam.baseeffects.df <- do.call(rbind, gam.baseeffects.df) #merge them into one 
  
  #Extract and combine interaction outputs
  gam.interactions.df <- lapply(gam.outputs.regionlist, '[[', "gam.covsmooth.interaction" ) #extract this df from each region's list
  gam.interactions.df <- do.call(rbind, gam.interactions.df) #merge them into one 
  
  gam.statistics.df <- list(gam.baseeffects.df, gam.interactions.df)
  names(gam.statistics.df) <- list("gam.baseeffects.df", "gam.interactions.df")
  saveRDS(gam.statistics.df, sprintf("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssociationsWithMyelin/%s", output.df.name))
}

colnames(SGIGR1_acw_electrodeatlas)[colnames(SGIGR1_acw_electrodeatlas) == "age.x"] <- "age"

R1.acw.interaction.gams(input.depth.df = SGIGR1_acw_electrodeatlas, acw.measure = "ACW_50", output.df.name = "R1ACW50_fix_depth_interaction.RDS")



## Plot the R1 ACW effects ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/statAssociationsWithMyelin/") 
files <- list.files(getwd()) 

#read in files and assign to variables
for(i in 1:length(files)){
  
  Rfilename <- gsub(".RDS", "", files[i])
  
  x <- readRDS(files[i]) 
  assign(Rfilename, x) 
}

R1ACW_50_superficial_fix_maineffect$gam.statistics.df
R1ACW_50_superficial_fix_maineffect$gam.statistics.df$depth <- 'superficial'


R1ACW_50_deep_fix_maineffect$gam.statistics.df
R1ACW_50_deep_fix_maineffect$gam.statistics.df$depth <- 'deep'


R1acw.plotdata <- rbind(R1ACW_50_superficial_fix_maineffect$gam.statistics.df , R1ACW_50_deep_fix_maineffect$gam.statistics.df)


R1acw.plotdata$orig_parcelname <- gsub("_R1", "", R1acw.plotdata$orig_parcelname)
R1acw.plotdata <- R1acw.plotdata %>% filter(orig_parcelname != "motor") #plot ones with significant depth interactions
R1acw.plotdata$orig_parcelname <- factor(R1acw.plotdata$orig_parcelname, ordered = T, levels = c("vlpfc", "dlpfc", "spfc")) #order for plotting

ggplot(R1acw.plotdata, aes (x = orig_parcelname, y = GAM.covsmooth.Fvalue, group = depth, color = interaction(depth))) +
  geom_point(size = 3.5, shape = 15) + 
  theme_classic() +
  scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6", "#969b98")) +
  labs(x="\n", y=sprintf("R1 - ACW_50 Statistic\n")) +
  theme(
    legend.position = "right", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


## R1 - acw interactions by cortical depth ----

#### Interaction effects

R1ACW50_fix_depth_interaction$gam.interactions.df %>% kbl() %>% kable_classic(full_width = F)

#### Marginal means 
vlpfc.interaction.model <- gam(vlpfc_R1 ~ s(age, k = 4, fx = F) + s(vlpfc_ACW_50, k = 3, fx = F) + s(vlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(vlpfc.interaction.model, terms = c("vlpfc_ACW_50[all]", "depth[all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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



dlpfc.interaction.model <- gam(dlpfc_R1 ~ s(age, k = 4, fx = F) + s(dlpfc_ACW_50, k = 3, fx = F) + s(dlpfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas)

interaction.plotdata <- ggpredict(dlpfc.interaction.model, terms = c("dlpfc_ACW_50 [all]", "depth [all]"), interval = "confidence")

plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


spfc.interaction.model <- gam(spfc_R1 ~ s(age, k = 4, fx = F) + s(spfc_ACW_50, k = 3, fx = F) + s(spfc_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(spfc.interaction.model, terms = c("spfc_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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


motor.interaction.model <- gam(motor_R1 ~ s(age, k = 4, fx = F) + s(motor_ACW_50, k = 3, fx = F) + s(motor_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), data = SGIGR1_acw_electrodeatlas) 

interaction.plotdata <- ggpredict(motor.interaction.model, terms = c("motor_ACW_50[all]", "depth[all]"), interval = "confidence")
plot(interaction.plotdata, show_residuals = T,  show_ci = F, colors = c("#012b0d", "#90b8a6"), line_size = 2, dot_size = 2, show_title = F) + 
  theme_classic() +
  labs(x="\nACW_50", y=sprintf("R1 (adjusted)\n")) +
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

