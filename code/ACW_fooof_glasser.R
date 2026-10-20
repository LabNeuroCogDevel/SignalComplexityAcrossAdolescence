
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
library(ggnewscale)

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()



lm.model.function <- function(input.df, region, dep_var, id_var){
  
  form <- as.formula(
    paste(dep_var, "~", id_var, "+ age + (1|lunaid)")
  )
  
  lm.model <- lmerTest::lmer(
    form,
    data = input.df %>%
      filter(orig_parcelname == region)
  )
  
  lm.stats <- data.frame(
    t(summary(lm.model)$coefficients[id_var, ])
  )
  return(lm.stats)
}

## ----------------------------------------------------------
## Load in and organize the glasser fooof files ----
## ----------------------------------------------------------

setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/fooof/")
regionOrder <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasserRegionOrder.csv')

# List all CSV files in the directory
csv_files <- list.files(pattern = "*AperiodicValues.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$channel <- sub(" @ .*", "", data$channel)
  
  
  avgData <- data %>%
    group_by(Subject, channel) %>%
    summarize(
      Exponent = mean(exponent, na.rm = TRUE),
      Offset = mean(offset, na.rm = TRUE),
      .groups = "drop"
    )

  
  combined_data <- rbind(combined_data, avgData)
}


write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/fooofMeasures_glasser.csv', row.names = F)



task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_taskRest_allregions_2.csv')

fooofAllChannels <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/fooofMeasures_glasser.csv')

fooofAllChannels <- fooofAllChannels %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", channel)) %>% select(-channel) %>% separate(Subject, c('lunaid','vdate'))


acw_fooof <- merge(fooofAllChannels, task_rest, by = c("lunaid", "vdate", "orig_parcelname"))




# --------------------------------------------------------------------
# EXPONENT ----
# --------------------------------------------------------------------



acw_fooof_rest <- acw_fooof %>% filter(epoch == 'Rest')

ggplot(acw_fooof_rest, aes(x = ACW_50)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "#27763d",
    color = "white"
  ) +
  geom_density(
    color = "#001a08",
    linewidth = 1
  ) +
  theme_classic()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = c("black")),
    axis.title.x = element_text(size = 12, color = c("black")),
    axis.title.y = element_text(size = 12, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) 


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_acw_hist.pdf", device = "pdf", dpi = 500, width = 2.3, height = 2.8)



ggplot(acw_fooof_rest, aes(x = Exponent)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "#27763d",
    color = "white"
  ) +
  geom_density(
    color = "#001a08",
    linewidth = 1
  ) +
  theme_classic()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = c("black")),
    axis.title.x = element_text(size = 12, color = c("black")),
    axis.title.y = element_text(size = 12, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) 


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_exponent_hist.pdf", device = "pdf", dpi = 500, width = 2.3, height = 2.8)



regions <- unique(acw_fooof_rest$orig_parcelname)
acw_fooof_rest$lunaid <- as.factor(acw_fooof_rest$lunaid)

# Initialize lists to store results
lmstats_fooof_list <- list()

# Loop through each combination of epoch and region
  for (loc in regions) {
    

    lm.stats <- lm.model.function(input.df = acw_fooof_rest, region = loc, dep_var = 'Exponent', id_var = 'ACW_50')
    
    # Create a data frame of stats and add epoch  info
    lm.stats  <-lm.stats  %>% 
      mutate(orig_parcelname = loc)
    
  
    # Store results
    lmstats_fooof_list[[loc]] <- lm.stats
   
  }



saveRDS(lmstats_fooof_list, file = "lmstats_list_fooof.rds")



lmstats_fooof_list <- readRDS("lmstats_list_fooof.rds")

lmstats_combined <- bind_rows(lmstats_fooof_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )


glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))


lmstats_combined <- lmstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)


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


lmstats_combined <- lmstats_combined %>%
  mutate(significant = Pr...t.. < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = Estimate)
  ) +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = Estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined,
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
  ) 


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/eastimateFOOOF_ACW_rest.pdf", device = "pdf", dpi = 500, width = 10, height = 5)



# --------------------------------------------------------------------
# Superior Temporal Sulcus ----
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_STSvp_ROI', 'L_STSvp_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Exponent = mean(Exponent, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Exponent, linetype=epoch, color = "#90b8a6"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#90b8a6")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Superior Temporal Sulcus (STSvp)') +
  xlab("ACW(s)") +
  ylab("Exponent")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_STS.pdf", device = "pdf", dpi = 500, width = 10, height = 10)



lm.model <- lmerTest::lmer("Exponent ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)


acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))



ggplot(acw_fooof_rest_avg, aes(x = Exponent)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "#27763d",
    color = "white"
  ) +
  geom_density(
    color = "#001a08",
    linewidth = 1
  ) +
  theme_classic()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = c("black")),
    axis.title.x = element_text(size = 12, color = c("black")),
    axis.title.y = element_text(size = 12, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Superior Temporal Sulcus (STSvp)')


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_exp_STS_hist.pdf", device = "pdf", dpi = 500, width = 2.3, height = 2.8)



# --------------------------------------------------------------------
# Visual
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_V4_ROI', 'L_V4_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Exponent = mean(Exponent, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Exponent, linetype=epoch, color = "#27763d"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#27763d")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Visual Cortex (V4)') +
  xlab("ACW(s)") +
  ylab("Exponent")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_V4.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Exponent ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))


ggplot(acw_fooof_rest_avg, aes(x = Exponent)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "#27763d",
    color = "white"
  ) +
  geom_density(
    color = "#001a08",
    linewidth = 1
  ) +
  theme_classic()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = c("black")),
    axis.title.x = element_text(size = 12, color = c("black")),
    axis.title.y = element_text(size = 12, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Visual Cortex (V4)')


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_exp_STS_hist.pdf", device = "pdf", dpi = 500, width = 2.3, height = 2.8)



# --------------------------------------------------------------------
# Motor
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Exponent = mean(Exponent, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Exponent, linetype=epoch, color = "#014421"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#014421")) + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Primary Motor Cortex (4)') +
  xlab("ACW(s)") +
  ylab("Exponent")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_motor.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Exponent ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))

# --------------------------------------------------------------------
# DLPFC ----
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_46_ROI', 'L_46_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Exponent = mean(Exponent, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Exponent, linetype=epoch, color = "#001a08"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3)+
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#001a08")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Dorsolateral Prefrontal Cortex (45)')+
  xlab("ACW(s)") +
  ylab("Exponent")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_dlpfc.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Exponent ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))

# --------------------------------------------------------------------
# OFFSET ----
# --------------------------------------------------------------------



acw_fooof_rest <- acw_fooof %>% filter(epoch == 'Rest')



ggplot(acw_fooof_rest, aes(x = Offset)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "#27763d",
    color = "white"
  ) +
  geom_density(
    color = "#001a08",
    linewidth = 1
  ) +
  theme_classic()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = c("black")),
    axis.title.x = element_text(size = 12, color = c("black")),
    axis.title.y = element_text(size = 12, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) 


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_offset_hist.pdf", device = "pdf", dpi = 500, width = 2.3, height = 2.8)


regions <- unique(acw_fooof_rest$orig_parcelname)
acw_fooof_rest$lunaid <- as.factor(acw_fooof_rest$lunaid)

# Initialize lists to store results
lmstats_fooof_list_offset <- list()

# Loop through each combination of epoch and region
for (loc in regions) {
  
  
  lm.stats <- lm.model.function(input.df = acw_fooof_rest, region = loc, dep_var = 'Offset', id_var = 'ACW_50')
  
  # Create a data frame of stats and add epoch  info
  lm.stats  <-lm.stats  %>% 
    mutate(orig_parcelname = loc)
  
  
  # Store results
  lmstats_fooof_list_offset[[loc]] <- lm.stats
  
}



saveRDS(lmstats_fooof_list_offset, file = "lmstats_list_fooof_offset.rds")

lmstats_fooof_list_offset <- readRDS("lmstats_list_fooof_offset.rds")

lmstats_combined_offset <- bind_rows(lmstats_fooof_list_offset, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )


glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))


lmstats_combined_offset <- lmstats_combined_offset %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)


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


lmstats_combined_offset <- lmstats_combined_offset %>%
  mutate(significant = Pr...t.. < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined_offset,
    atlas = glasser,
    side = "lateral",
    aes(fill = Estimate)
  ) +
  geom_brain(
    data = lmstats_combined_offset,
    atlas = glasser,
    side = "lateral",
    aes(fill = Estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined_offset,
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
  ) 


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/eastimateFOOOF_ACW_rest_offset.pdf", device = "pdf", dpi = 500, width = 10, height = 5)



# --------------------------------------------------------------------
# Superior Temporal Sulcus ----
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_STSvp_ROI', 'L_STSvp_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Offset = mean(Offset, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Offset, linetype=epoch, color = "#90b8a6"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#90b8a6")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Superior Temporal Sulcus (STSvp)') +
  xlab("ACW(s)") +
  ylab("Offset")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_STS_offset.pdf", device = "pdf", dpi = 500, width = 10, height = 10)


lm.model <- lmerTest::lmer("Offset ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)



acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))

# --------------------------------------------------------------------
# Visual ----
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_V4_ROI', 'L_V4_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Offset = mean(Offset, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg , aes(x = ACW_50, y = Offset, linetype=epoch, color = "#27763d"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#27763d")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Visual Cortex (V4)') +
  xlab("ACW(s)") +
  ylab("Offset")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_V4_offset.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Offset ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))

# --------------------------------------------------------------------
# Motor ----
# --------------------------------------------------------------------
acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Offset = mean(Offset, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Offset, linetype=epoch, color = "#014421"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#014421")) + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Primary Motor Cortex (4)') +
  xlab("ACW(s)") +
  ylab("Offset")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_motor_offset.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Offset ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))

# --------------------------------------------------------------------
# DLPFC ----
# --------------------------------------------------------------------

acw_fooof_rest_avg <- acw_fooof_rest %>%
  filter(orig_parcelname %in% c('R_46_ROI', 'L_46_ROI')) %>%
  group_by(lunaid, vdate, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    Offset = mean(Offset, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = acw_fooof_rest_avg, aes(x = ACW_50, y = Offset, linetype=epoch, color = "#001a08"))) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 1.5) +
  geom_point(alpha=0.3)+
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  scale_color_manual(values = c("#001a08")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + ggtitle('Dorsolateral Prefrontal Cortex (45)')+
  xlab("ACW(s)") +
  ylab("Offset")

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure2/FOOOF_ACW_dlpfc_offset.pdf", device = "pdf", dpi = 500, width = 10, height = 10)

lm.model <- lmerTest::lmer("Offset ~ ACW_50 + age + (1|lunaid)", data = acw_fooof_rest_avg)

summary(lm.model)

acw_fooof_rest_avg %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))
