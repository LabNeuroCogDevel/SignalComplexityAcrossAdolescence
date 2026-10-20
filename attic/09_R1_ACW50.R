library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(kableExtra)
conflicted::conflicts_prefer(dplyr::filter)

#EEG electrode atlas cortical areas list

electrodes <- c("vlpfc", "dlpfc", "spfc", "motor")

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


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# Task, Rest, vs R1 Frontal, Parietal, Occipital ----

electrodes <- c("Frontal", "Parietal", "Occipital")

SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/superDeep_frontalParietalOccipital_R1.RDS")

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth$subject_id <- as.factor(depth$subject_id) #factor for gam modeling
  return(depth)
})

SGIGR1_acw_electrodeatlas <- do.call(rbind, SGIGmyelin.EEGatlas.7T)
SGIGR1_acw_electrodeatlas$depth <- factor(sub("\\..*$", "", row.names(SGIGR1_acw_electrodeatlas)), levels = c("deep", "superficial"), ordered = F) 

SGIGR1_acw_electrodeatlas <- SGIGR1_acw_electrodeatlas %>% group_by(depth) %>%
  mutate(across(c("Frontal_R1","Parietal_R1", "Occipital_R1"), naoutlier))


task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/taskRest.csv')

task_rest_wide <- task_rest %>%
  select(-X) %>%  # Remove the unnecessary column "X"
  pivot_wider(names_from = region, values_from = ACW_50, names_glue = "{region}_ACW_50")

R1_ACW_FPO_task <- merge(SGIGR1_acw_electrodeatlas, task_rest_wide %>% filter(epoch == 'Task'), by = c('lunaid', 'visitno', 'sex'))

R1_ACW_FPO_rest <- merge(SGIGR1_acw_electrodeatlas, task_rest_wide %>% filter(epoch == 'restEyesOpen'), by = c('lunaid', 'visitno', 'sex'))

colnames(R1_ACW_FPO_task)[colnames(R1_ACW_FPO_task) == "age.x"] <- "age"
colnames(R1_ACW_FPO_rest)[colnames(R1_ACW_FPO_rest) == "age.x"] <- "age"


## TASK ----
R1.acw.maineffect.gams(input.depth.df = R1_ACW_FPO_task %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_FPO_task_maineffect.RDS")
R1.acw.maineffect.gams(input.depth.df = R1_ACW_FPO_task %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_FPO_task_maineffect.RDS")

#### Marginal means 
frontal.interaction.model <- gam(Frontal_R1 ~ s(age, k = 4, fx = F)  + s(Frontal_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                               data = R1_ACW_FPO_task)

summary(frontal.interaction.model)
p.adjust(0.04, method = 'bonferroni', n=3)

interaction.plotdata <- ggpredict(frontal.interaction.model, terms = c("Frontal_ACW_50[all]", "depth[all]"), interval = "confidence")

lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_task, aes(x = Frontal_ACW_50, y = Frontal_R1, shape=depth, color = '#012b0d'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#012b0d'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#012b0d"), guide='none') + 
          scale_fill_manual(values = c("#012b0d"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Frontal_myelin_acw_task.pdf", device = "pdf", dpi = 500, width = 6, height = 6)



Parietal.interaction.model <- gam(Parietal_R1 ~ s(age, k = 4, fx = F) + s(Parietal_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                                  data = R1_ACW_FPO_task)
summary(Parietal.interaction.model)

interaction.plotdata <- ggpredict(Parietal.interaction.model, terms = c("Parietal_ACW_50 [all]", "depth [all]"), interval = "confidence")

lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_task, aes(x = Parietal_ACW_50, y = Parietal_R1, shape=depth, color = '#27763d'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#27763d'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#27763d"), guide='none') + 
          scale_fill_manual(values = c("#27763d"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Parietal_myelin_acw_task.pdf", device = "pdf", dpi = 500, width = 6, height = 6)



Occipital.interaction.model <- gam(Occipital_R1 ~ s(age, k = 4, fx = F) + s(Occipital_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                                   data = R1_ACW_FPO_task) 
summary(Occipital.interaction.model)

interaction.plotdata <- ggpredict(Occipital.interaction.model, terms = c("Occipital_ACW_50[all]", "depth[all]"), interval = "confidence")


lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_task, aes(x = Occipital_ACW_50, y = Occipital_R1, shape=depth, color = '#90b8a6'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#90b8a6'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#90b8a6"), guide='none') + 
          scale_fill_manual(values = c("#90b8a6"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Occipital_myelin_acw_task.pdf", device = "pdf", dpi = 500, width = 6, height = 6)

## REST ----
R1.acw.maineffect.gams(input.depth.df = R1_ACW_FPO_rest %>% filter(depth == 'superficial'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_superficial_FPO_rest_maineffect.RDS")
R1.acw.maineffect.gams(input.depth.df = R1_ACW_FPO_rest %>% filter(depth == 'deep'), acw.measure = "ACW_50", output.df.name = "R1ACW_50_deep_FPO_rest_maineffect.RDS")

#### Marginal means 
frontal.interaction.model <- gam(Frontal_R1 ~ s(age, k = 4, fx = F)  + s(Frontal_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                                 data = R1_ACW_FPO_rest)

summary(frontal.interaction.model)
p.adjust(0.0177, method = 'bonferroni', n=3)

interaction.plotdata <- ggpredict(frontal.interaction.model, terms = c("Frontal_ACW_50[all]", "depth[all]"), interval = "confidence")


lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_rest, aes(x = Frontal_ACW_50, y = Frontal_R1, shape=depth, color = '#012b0d'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#012b0d'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#012b0d"), guide='none') + 
          scale_fill_manual(values = c("#012b0d"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Frontal_myelin_acw_rest.pdf", device = "pdf", dpi = 500, width = 6, height = 6)



Parietal.interaction.model <- gam(Parietal_R1 ~ s(age, k = 4, fx = F) + s(Parietal_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                                  data = R1_ACW_FPO_rest)
summary(Parietal.interaction.model)

interaction.plotdata <- ggpredict(Parietal.interaction.model, terms = c("Parietal_ACW_50 [all]", "depth [all]"), interval = "confidence")

lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_rest, aes(x = Parietal_ACW_50, y = Parietal_R1, shape=depth, color = '#27763d'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#27763d'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#27763d"), guide='none') + 
          scale_fill_manual(values = c("#27763d"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Parietal_myelin_acw_rest.pdf", device = "pdf", dpi = 500, width = 6, height = 6)



Occipital.interaction.model <- gam(Occipital_R1 ~ s(age, k = 4, fx = F) + s(Occipital_ACW_50, by = depth, k = 3) + s(subject_id, bs="re") + depth, method = c("REML"), 
                                   data = R1_ACW_FPO_rest) 
summary(Occipital.interaction.model)

interaction.plotdata <- ggpredict(Occipital.interaction.model, terms = c("Occipital_ACW_50[all]", "depth[all]"), interval = "confidence")

lunaize(ggplot( ) + 
          geom_point(data = R1_ACW_FPO_rest, aes(x = Occipital_ACW_50, y = Occipital_R1, shape=depth, color = '#90b8a6'), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, linetype=group, color='#90b8a6'), data=interaction.plotdata, size = 2)+ 
          theme_minimal()+ 
          scale_color_manual(values = c("#90b8a6"), guide='none') + 
          scale_fill_manual(values = c("#90b8a6"), guide='none') + 
          labs(x = "ACW", y = "R1 (adjusted)")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Occipital_myelin_acw_rest.pdf", device = "pdf", dpi = 500, width = 6, height = 6)


