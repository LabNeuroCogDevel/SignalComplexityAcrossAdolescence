
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
library(extrafont)
extrafont::loadfonts()
library(conflicted)
library(eegkit)
library(ggseg)
library(ggsegGlasser)

dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#ffffff"   # white highlight
)


task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_taskRest_allregions_2.csv')

meanINTperPerson <- task_rest %>%
  group_by(lunaid, visitno, age, sex, epoch) %>%
  summarise(meanINT = mean(ACW_50, na.rm=T)) %>%
  ungroup()


lunaize(ggplot(data = meanINTperPerson, aes(x = age, y = meanINT, color = epoch))) + 
  geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2) +
  geom_smooth(aes(group=epoch), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) + 
  scale_color_manual(values = c("#001a08", "#90b8a6")) + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) +
  ggtitle('MeanINT')




task_rest_centered <- task_rest %>%
  group_by(lunaid, visitno, age, sex, epoch) %>%
  mutate(
    subject_mean = mean(ACW_50, na.rm = TRUE),
    ACW_50_centered = ACW_50 - subject_mean
  ) %>%
  ungroup()



lunaize(ggplot(data = task_rest_centered, aes(x = age, y = ACW_50_centered, color=epoch))) + 
  geom_smooth(aes(group=interaction(lobe, epoch)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.03, linewidth = 1.5) + 
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) +
  ggtitle('mean normalized') +
  facet_wrap(~lobe)





## Principal developmental gradient
epochs <- unique(task_rest_centered$epoch)
regions <- unique(task_rest_centered$orig_parcelname)
task_rest_centered$lunaid <- as.factor(task_rest_centered$lunaid)
task_rest_centered$epoch <- as.factor(task_rest_centered$epoch)

# Initialize lists to store results
gamstats_glasser_list <- list()
gamderivatives_glasser_list <- list()
gamsmoothestimates_glasser_list <- list()
gamfittedValues_glasser_list <- list()

# Loop through each combination of epoch and region
for (task in epochs) {
  for (loc in regions) {
    
    # Apply the growthrate function for each combination
    allstats <- growthrate(task_rest_centered %>% dplyr::filter(epoch == task & orig_parcelname == loc), 'ACW_50', 'age', loc)
    
    # Create a data frame of stats and add epoch  info
    gamstats <- allstats$gam.statistics %>% 
      mutate(epoch = task) %>%
      mutate(orig_parcelname = loc)
    
    # Create a data frame of derivatives and add epoch info
    gamderivatives <- allstats$gam.derivatives %>% 
      mutate(epoch = task) %>%
      mutate(orig_parcelname = loc)
    
    gamsmoothestimates <- allstats$gam.smoothestimates %>% 
      mutate(epoch = task) %>%
      mutate(orig_parcelname = loc)
    
    gamfittedValues <- allstats$gam.fittedvalues %>% 
      mutate(epoch = task) %>%
      mutate(orig_parcelname = loc)
    
    # Store results
    gamstats_glasser_list[[paste(task,loc, sep = "_")]] <- gamstats
    gamderivatives_glasser_list[[paste(task, loc, sep = "_")]] <- gamderivatives
    gamsmoothestimates_glasser_list[[paste(task,loc, sep = "_")]] <- gamsmoothestimates
    gamfittedValues_glasser_list[[paste(task,loc, sep = "_")]] <- gamfittedValues
    
    
  }
}

saveRDS(gamstats_glasser_list, file = "gamstats_list_allregions.rds")
saveRDS(gamderivatives_glasser_list, file = "gamderivatives_list_allregions.rds")
saveRDS(gamsmoothestimates_glasser_list, file = "gamsmoothestimates_list_allregions.rds")
saveRDS(gamfittedValues_glasser_list, file = "gamfittedValues_list_allregions.rds")


# ----------------------------------
## PCA on glasser trajectories
# ----------------------------------
gamstats_meanCentered_list <- readRDS("gamstats_list_allregions.rds")
gamderivatives_meanCentered_list <- readRDS("gamderivatives_list_allregions.rds")

gamstats_combined_meanCentered <- bind_rows(gamstats_meanCentered_list, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

gamderivatives_combined_meanCentered <- bind_rows(gamderivatives_meanCentered_list, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

gamstats_combined_meanCentered <- merge(gamstats_combined_meanCentered, glasser_lobes, by = c('orig_parcelname'))



dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#ffffff"   # white highlight
)

ggplot() +
  theme_void() +
  geom_brain(
    data = gamstats_combined_meanCentered,
    atlas = glasser,
    side = "lateral",
    aes(fill = mean.derivative)
  ) +
  scale_fill_gradientn(
    colours =  dramatic_palette) + facet_wrap(~epoch)

ggplot() +
  theme_void() +
  geom_brain(
    data = gamstats_combined %>% dplyr::filter(epoch == "Task"),
    atlas = glasser,
    side = "lateral",
    aes(fill = mean.derivative)
  ) +
  scale_fill_gradientn(
    colours =  dramatic_palette)


gam.smoothestimates.glasser <- bind_rows(gamsmoothestimates_glasser_list, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

gam.smoothestimates.glasser <- gam.smoothestimates.glasser %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

  
gam.smoothestimates.glasser.rest <- gam.smoothestimates.glasser %>% dplyr::filter(epoch == 'Rest') %>% select(label, age, .estimate)
gam.smoothestimates.rest.long <- gam.smoothestimates.glasser.rest %>% pivot_wider(names_from="label",values_from=".estimate")
gam.smoothestimates.rest.long <- gam.smoothestimates.rest.long %>% select(-age)  #ages by regions matrix for PCA

smooths.pca.rest <- prcomp(gam.smoothestimates.rest.long) #PCA

summary(smooths.pca.rest)$importance[2,1]


PC1 <- as.data.frame(smooths.pca.rest$rotation[,1]) #first principal component loadings
PC1$label <- row.names(PC1)
colnames(PC1) <- c("PC1","label")
PC1$ranked <- rank(PC1$PC1)

trajectory_df <- data.frame(
  age = sort(unique(gam.smoothestimates.glasser.rest$age)),
  PC1 = smooths.pca.rest$x[,1],
  PC2 = smooths.pca.rest$x[,2]
)

ggplot(trajectory_df, aes(age, PC1)) +
  geom_line(linewidth = 1.2) +
  theme_classic()



# Plot on cortex
ggplot() +
  theme_void() +
  geom_brain(
    data = PC1,
    atlas = glasser,
    side = "lateral",
    aes(fill = PC1)
  ) +
  scale_fill_gradient(
    low = "#d8efe2",   # light
    high = "#012b0d",  # dark
    na.value = "grey90"
  )





gam.smoothestimates.glasser.task <- gam.smoothestimates.glasser %>% dplyr::filter(epoch == 'Task') %>% select(label, age, .estimate)
gam.smoothestimates.task.long <- gam.smoothestimates.glasser.task %>% pivot_wider(names_from="label",values_from=".estimate")
gam.smoothestimates.task.long <- gam.smoothestimates.task.long %>% select(-age)  #ages by regions matrix for PCA

smooths.pca.task <- prcomp(gam.smoothestimates.task.long) #PCA

summary(smooths.pca.task)$importance[2,1]


PC1_task <- as.data.frame(smooths.pca.task$rotation[,1]) #first principal component loadings
PC1_task$label <- row.names(PC1_task)
colnames(PC1_task) <- c("PC1","label")
PC1_task$ranked <- rank(PC1_task$PC1)


trajectory_df_task <- data.frame(
  age = sort(unique(gam.smoothestimates.glasser.task$age)),
  PC1 = smooths.pca.rest$x[,1],
  PC2 = smooths.pca.rest$x[,2]
)


ggplot(trajectory_df_task, aes(age, PC1)) +
  geom_line(linewidth = 1.2) +
  theme_classic()



# Plot on cortex
ggplot() +
  theme_void() +
  geom_brain(
    data = PC1_task,
    atlas = glasser,
    side = "lateral",
    aes(fill = PC1)
  ) +
  scale_fill_gradient(
    low = "#d8efe2",   # light
    high = "#012b0d",  # dark
    na.value = "grey90"
  )


cor(PC1_task$PC1, PC1$PC1)
