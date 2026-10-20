
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
library(eegkit)
library(ggseg)
library(ggsegGlasser)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::mutate)
library(sf)

# ----------------------------------------
# load eeg data 
# ----------------------------------------

task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_taskRest_allregions_2.csv')
names(task_rest)[names(task_rest) == "vdate"] <- "rest.date"

# ----------------------------------------
# load and prep myelin data 
# ----------------------------------------

myelin.compartments.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/compartmentsR1_glasseratlas_finalsample.RDS")
glasser.snr.exclude <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/SNR/glasser_SNR_exclusion.csv")

myelin.compartments.7T.long <- lapply(myelin.compartments.7T, function(depth){
  cols_to_pivot <- names(depth)[grepl("ROI", names(depth))]
  depth.long <- depth %>% pivot_longer(cols = cols_to_pivot, names_to = "orig_parcelname", values_to = "R1")
  depth.long <- depth.long[!(depth.long$orig_parcelname %in% glasser.snr.exclude$orig_parcelname),]
  depth.long <- depth.long %>% mutate(subses = sprintf("%s_%s", subject_id, session_id))
  return(depth.long)
})

myelin.compartments.7T.long.superficial <- myelin.compartments.7T.long$superficial

myelin.compartments.7T.long.middle <- myelin.compartments.7T.long$middle

myelin.compartments.7T.long.deep <- myelin.compartments.7T.long$deep


# ----------------------------------------
# merge myelin and eeg by depth
# ----------------------------------------


myelin.eeg.superficial <- merge(myelin.compartments.7T.long.superficial, task_rest %>% select(-age), by = c('lunaid', 'rest.date', 'orig_parcelname', 'visitno', 'sex'))
myelin.eeg.middle <- merge(myelin.compartments.7T.long.middle, task_rest %>% select(-age), by = c('lunaid', 'rest.date', 'orig_parcelname', 'visitno', 'sex'))
myelin.eeg.deep <- merge(myelin.compartments.7T.long.deep, task_rest %>% select(-age), by = c('lunaid', 'rest.date', 'orig_parcelname', 'visitno', 'sex'))



# ----------------------------------------
# eeg myelin gams - Superficial 
# ----------------------------------------


epochs <- unique(myelin.eeg.superficial$epoch)
regions <- unique(myelin.eeg.superficial$orig_parcelname)
myelin.eeg.superficial$lunaid <- as.factor(myelin.eeg.superficial$lunaid)
myelin.eeg.superficial$epoch <- as.factor(myelin.eeg.superficial$epoch)

# Initialize lists to store results
gamstats_glasser_list <- list()
gamderivatives_glasser_list <- list()
gamsmoothestimates_glasser_list <- list()
gamfittedValues_glasser_list <- list()


# Loop through each combination of epoch and region
for (task in epochs) {
  for (loc in regions) {
    
    # Apply the growthrate function for each combination
    allstats <- gam.statistics.smooths(input.df = myelin.eeg.superficial %>% dplyr::filter(epoch == task), region = loc, dep_var = 'ACW_50', 
                                       smooth_var = "age", id_var = "lunaid", covariates = "R1", 
                                       random_intercepts = TRUE, random_slopes = FALSE, knots = 4, set_fx = FALSE)
    
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

saveRDS(gamstats_glasser_list, file = "gamstats_list_EEGmyelin.rds")
saveRDS(gamderivatives_glasser_list, file = "gamderivatives_list_EEGmyelin.rds")
saveRDS(gamsmoothestimates_glasser_list, file = "gamsmoothestimates_list_EEGmyelin.rds")
saveRDS(gamfittedValues_glasser_list, file = "gamfittedValues_list_EEGmyelin.rds")



gamstats_glasser_list <- readRDS("gamstats_list_EEGmyelin.rds")
gamderivatives_glasser_list <- readRDS("gamderivatives_list_EEGmyelin.rds")

gamstats_combined <- bind_rows(gamstats_glasser_list, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

gamderivatives_combined <- bind_rows(gamderivatives_glasser_list, .id = "epoch_region") %>%
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

gamderivatives_combined <- gamderivatives_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)



gamstats_combined <- gamstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

glasser.regions <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/HCPMMP_glasseratlas/glasser360_regionlist.csv")

glasser.plotting <- glasser.regions
glasser.plotting$cortex <- "cortex"
glasser.plotting$cortex[glasser.plotting$orig_parcelname %in% glasser.snr.exclude$orig_parcelname] <- "exclude"
glasser.plotting <- glasser.plotting %>% select(label, cortex)
glasser.plotting <- rbind(glasser.plotting, data.frame(label = "rh_???", cortex = "exclude"))
glasser.plotting <- rbind(glasser.plotting, data.frame(label = "lh_???", cortex = "exclude"))
glasser.plotting <- glasser.plotting %>% filter(label != "lh_L_10pp")



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


gamstats_combined <- gamstats_combined %>%
  mutate(significant = GAM.smooth.pvalue < 0.05)



ggplot() +
  theme_void() +
  geom_brain(
    data = gamstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = GAM.smooth.Fvalue)
  ) +
  geom_brain(
    data = gamstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = GAM.smooth.Fvalue),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  (dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = gamstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) + 
  new_scale_fill() +
  geom_brain(data = glasser.plotting, 
             atlas = glasser, 
             side = "lateral",
             mapping = aes(fill = cortex), 
             colour=I(alpha("gray20", 0)), size=I(0)) +
  scale_fill_manual(values = c(alpha("white", 0), "gray75")) + facet_wrap(~epoch)


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/INT_superficialMyelin.pdf", device = "pdf", dpi = 500, width = 18, height = 5)





# --------------------------------------------------------------------
# Superior Temporal Sulcus
# --------------------------------------------------------------------
myelin.eeg.superficial_avg <- myelin.eeg.superficial %>%
  filter(orig_parcelname %in% c('R_STSvp_ROI', 'L_STSvp_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.superficial_avg, aes(y = ACW_50, x = R1, linetype=epoch, color = "#90b8a6"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#90b8a6")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#90b8a6'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('STS (STSvp)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_STS.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)


gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.superficial_avg, random=list(lunaid=~1))
summary(gam.model$gam)



# --------------------------------------------------------------------
## Visual
# --------------------------------------------------------------------
myelin.eeg.superficial_avg <- myelin.eeg.superficial %>%
  filter(orig_parcelname %in% c('R_V4_ROI', 'L_V4_ROI')) %>%
  group_by(lunaid, rest.date, epoch,age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.superficial_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#27763d"))) +
  geom_smooth(aes(group=epoch), method='gam', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#27763d")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#27763d'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Visual Cortex (V4)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_V1.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)

gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.superficial_avg, random=list(lunaid=~1))
summary(gam.model$gam)




# --------------------------------------------------------------------
# Motor
# --------------------------------------------------------------------
myelin.eeg.superficial_avg <- myelin.eeg.superficial %>%
  filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.superficial_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#014421"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#014421")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#014421'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Motor Cortex (4)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_motor.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)


gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.superficial_avg, random=list(lunaid=~1))
summary(gam.model$gam)



# --------------------------------------------------------------------
# dlpfc
# --------------------------------------------------------------------
myelin.eeg.superficial_avg <- myelin.eeg.superficial %>%
  filter(orig_parcelname %in% c('R_46_ROI', 'L_46_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.superficial_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#001a08"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#001a08")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#001a08'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Motor Cortex (4)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_dlpfc.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)

gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.superficial_avg, random=list(lunaid=~1))
summary(gam.model$gam)

# ----------------------------------------
# eeg myelin gams - Deep 
# ----------------------------------------


epochs <- unique(myelin.eeg.deep$epoch)
regions <- unique(myelin.eeg.deep$orig_parcelname)
myelin.eeg.deep$lunaid <- as.factor(myelin.eeg.deep$lunaid)
myelin.eeg.deep$epoch <- as.factor(myelin.eeg.deep$epoch)

# Initialize lists to store results
gamstats_glasser_list_deep <- list()
gamderivatives_glasser_list_deep <- list()
gamsmoothestimates_glasser_list_deep <- list()
gamfittedValues_glasser_list_deep <- list()


# Loop through each combination of epoch and region
for (task in epochs) {
  for (loc in regions) {
    
    # Apply the growthrate function for each combination
    allstats <- gam.statistics.smooths(input.df = myelin.eeg.deep %>% dplyr::filter(epoch == task), region = loc, dep_var = 'ACW_50', 
                                       smooth_var = "age", id_var = "lunaid", covariates = "R1", 
                                       random_intercepts = TRUE, random_slopes = FALSE, knots = 4, set_fx = FALSE)
    
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
    gamstats_glasser_list_deep[[paste(task,loc, sep = "_")]] <- gamstats
    gamderivatives_glasser_list_deep[[paste(task, loc, sep = "_")]] <- gamderivatives
    gamsmoothestimates_glasser_list_deep[[paste(task,loc, sep = "_")]] <- gamsmoothestimates
    gamfittedValues_glasser_list_deep[[paste(task,loc, sep = "_")]] <- gamfittedValues
    
    
  }
}

saveRDS(gamstats_glasser_list_deep, file = "gamstats_list_EEGmyelin_deep.rds")
saveRDS(gamderivatives_glasser_list_deep, file = "gamderivatives_list_EEGmyelin_deep.rds")
saveRDS(gamsmoothestimates_glasser_list_deep, file = "gamsmoothestimates_list_EEGmyelin_deep.rds")
saveRDS(gamfittedValues_glasser_list_deep, file = "gamfittedValues_list_EEGmyelin_deep.rds")



gamstats_glasser_list_deep <- readRDS("gamstats_list_EEGmyelin_deep.rds")
gamderivatives_glasser_list_deep <- readRDS("gamderivatives_list_EEGmyelin_deep.rds")

gamstats_combined_deep <- bind_rows(gamstats_glasser_list_deep, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

gamderivatives_combined_deep <- bind_rows(gamderivatives_glasser_list_deep, .id = "epoch_region") %>%
  separate(
    epoch_region,
    into = c("epoch", "orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )
gamderivatives_combined_deep <- gamderivatives_combined_deep %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)


glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

gamstats_combined_deep <- gamstats_combined_deep %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)


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
gamstats_combined_deep <- gamstats_combined_deep %>%
  mutate(significant = GAM.smooth.pvalue < 0.05)



ggplot() +
  theme_void() +
  geom_brain(
    data = gamstats_combined_deep,
    atlas = glasser,
    side = "lateral",
    aes(fill = GAM.smooth.Fvalue)
  ) +
  geom_brain(
    data = gamstats_combined_deep,
    atlas = glasser,
    side = "lateral",
    aes(fill = GAM.smooth.Fvalue),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  (dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = gamstats_combined_deep,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) + 
  new_scale_fill() +
  geom_brain(data = glasser.plotting, 
             atlas = glasser, 
             side = "lateral",
             mapping = aes(fill = cortex), 
             colour=I(alpha("gray20", 0)), size=I(0)) +
  scale_fill_manual(values = c(alpha("white", 0), "gray75")) + facet_wrap(~epoch)


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/INT_deepMyelin.pdf", device = "pdf", dpi = 500, width = 18, height = 5)





# --------------------------------------------------------------------
# Superior Temporal Sulcus
# --------------------------------------------------------------------
myelin.eeg.deep_avg <- myelin.eeg.deep %>%
  filter(orig_parcelname %in% c('R_STSvp_ROI', 'L_STSvp_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.deep_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#90b8a6"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#90b8a6")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#90b8a6'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('STS (STSvp)') +
  ylab("ACW(s)") +
  xlab("R1")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_STS_deep.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)

gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.deep_avg, random=list(lunaid=~1))
summary(gam.model$gam)

# --------------------------------------------------------------------
## Visual
# --------------------------------------------------------------------
myelin.eeg.deep_avg <- myelin.eeg.deep %>%
  filter(orig_parcelname %in% c('R_V4_ROI', 'L_V4_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.deep_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#27763d"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#27763d")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#27763d'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Visual Cortex (V4)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_V1_deep.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)

gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.deep_avg, random=list(lunaid=~1))
summary(gam.model$gam)



# --------------------------------------------------------------------
# Motor
# --------------------------------------------------------------------
myelin.eeg.deep_avg <- myelin.eeg.deep %>%
  filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.deep_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#014421"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#014421")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#014421'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Motor Cortex (4)') +
  ylab("ACW(s)") +
  xlab("R1")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_motor_deep.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)


gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.deep_avg, random=list(lunaid=~1))
summary(gam.model$gam)



# --------------------------------------------------------------------
# dlpfc
# --------------------------------------------------------------------
myelin.eeg.deep_avg <- myelin.eeg.deep %>%
  filter(orig_parcelname %in% c('R_46_ROI', 'L_46_ROI')) %>%
  group_by(lunaid, rest.date, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    R1 = mean(R1, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = myelin.eeg.deep_avg, aes(x = ACW_50, y = R1, linetype=epoch, color = "#001a08"))) +
  geom_smooth(aes(group=epoch), method='lm', alpha=0.3, linewidth = 0.3) +
  geom_point(alpha=0.2 , aes(shape = epoch), size = 0.2) +
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth = 0.1) +
  scale_color_manual(values = c("#001a08")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 8, color = '#001a08'),
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .15), 
    axis.ticks = element_line(linewidth = .15),
    strip.background = element_rect('white')) + ggtitle('Dorsolateral PFC (46)') +
  ylab("ACW(s)") +
  xlab("R1")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure3/R1_ACW_dlpfc_deep.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 1.8, height = 2)


gam.model <-  mgcv::gamm(R1 ~ s(age, k = 3)  + epoch, data = myelin.eeg.deep_avg, random=list(lunaid=~1))
summary(gam.model$gam)
