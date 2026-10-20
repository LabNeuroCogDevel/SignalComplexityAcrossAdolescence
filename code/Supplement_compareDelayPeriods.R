
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

conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# ACW ----

allSub_delayPeriods <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_delayPeriods.csv')

# Count the number of unique delay lengths per subject
subject_delay_counts <- allSub_delayPeriods %>%
  group_by(SubjectID) %>%
  summarise(n_delays = n_distinct(delayLength))

# Find subjects missing all 3 delay lengths
subjects_missing_all <- subject_delay_counts %>%
  filter(n_delays != 3)

subjects_missing_all


allSub_delayPeriods_outlier <- allSub_delayPeriods %>% group_by(SubjectID, delayLength, lobe) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) 



lobe_avg <- allSub_delayPeriods_outlier %>%
  group_by(SubjectID, delayLength, lobe, age, sex, visitno) %>%  # keep these in group
  summarise(
    ACW_50_mean = mean(ACW_50, na.rm = TRUE),
    .groups = "drop"
  )%>% separate(SubjectID, c('lunaID','vdate'))



lunaize(
  ggplot(data = lobe_avg, 
         aes(x = delayLength, y = ACW_50_mean, fill = interaction(delayLength), color=delayLength)) +
    geom_violin(alpha = 1, position = position_dodge(width = 0.75), color = "black") +
    scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6")) + 
    scale_fill_manual(values = c("#012b0d", "#27763d", "#90b8a6")) + 
    stat_summary(
      fun = "mean", 
      geom = "point", 
      aes(group = interaction(delayLength)),  # Ensure grouping by ageGroup
      position = position_dodge(width = 0.75), 
      shape = 23, size = 3, fill = "white", color = "black"
    )+
    labs(x = "", y = "ACW", fill = "") +
    theme_minimal() + facet_wrap(~lobe, ncol=4) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')))


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Supp/delayLengths.pdf", device = "pdf", dpi = 500, width = 10, height = 4)



# Make delayLength an ordered factor
lobe_avg$delayLength_order <- factor(lobe_avg$delayLength, 
                               levels = c(6, 8, 10),  # specify the order
                               ordered = TRUE)

gam.model <-  gam(ACW_50_mean ~ s(age, k = 3) + s(age, k=3, by=delayLength_order) + delayLength + sex, data = lobe_avg %>% filter(lobe == 'F'), random=list(lunaID=~1))
summary(gam.model)

draw(gam.model)

lunaize(
  ggplot(data = lobe_avg, 
         aes(x = age, y = ACW_50_mean, color=delayLength)) + geom_point() +
    geom_smooth(aes(group=delayLength), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
    labs(x = "", y = "ACW", fill = "") +
    theme_minimal()
) + 
  theme(
    legend.position = "left",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~lobe)




#avg the delay periods together

delay_avg <- allSub_delayPeriods_outlier %>%
  group_by(SubjectID, orig_parcelname, lobe, age, sex, visitno) %>%  # keep these in group
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    .groups = "drop"
  )%>% separate(SubjectID, c('lunaID','vdate'))

write.csv(delay_avg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam.csv', row.names = F)



# Variability ----

allSub_delayPeriods <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_delayPeriods_variability.csv')

# Count the number of unique delay lengths per subject
subject_delay_counts <- allSub_delayPeriods %>%
  group_by(SubjectID) %>%
  summarise(n_delays = n_distinct(delayLength))

# Find subjects missing all 3 delay lengths
subjects_missing_all <- subject_delay_counts %>%
  filter(n_delays != 3)

subjects_missing_all


allSub_delayPeriods_outlier <- allSub_delayPeriods %>% group_by(SubjectID, delayLength, lobe) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) 



lobe_avg <- allSub_delayPeriods_outlier %>%
  group_by(SubjectID, delayLength, lobe, age, sex, visitno) %>%  # keep these in group
  summarise(
    ACW_50_mean = mean(ACW_50, na.rm = TRUE),
    .groups = "drop"
  )%>% separate(SubjectID, c('lunaID','vdate'))



lunaize(
  ggplot(data = lobe_avg, 
         aes(x = delayLength, y = ACW_50_mean, fill = interaction(delayLength), color=delayLength)) +
    geom_violin(alpha = 1, position = position_dodge(width = 0.75), color = "black") +
    scale_color_manual(values = c("#012b0d", "#27763d", "#90b8a6")) + 
    scale_fill_manual(values = c("#012b0d", "#27763d", "#90b8a6")) + 
    stat_summary(
      fun = "mean", 
      geom = "point", 
      aes(group = interaction(delayLength)),  # Ensure grouping by ageGroup
      position = position_dodge(width = 0.75), 
      shape = 23, size = 3, fill = "white", color = "black"
    )+
    labs(x = "", y = "ACW", fill = "") +
    theme_minimal() + facet_wrap(~lobe) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 20, color = c("black")),
      axis.title.x = element_text(size = 20, color = c("black")),
      axis.title.y = element_text(size = 20, color = c("black")),
      axis.line = element_line(linewidth = .2), 
      axis.ticks = element_line(linewidth = .2),
      strip.background = element_rect('white')))


# Make delayLength an ordered factor
lobe_avg$delayLength_order <- factor(lobe_avg$delayLength, 
                                     levels = c(6, 8, 10),  # specify the order
                                     ordered = TRUE)

gam.model <-  gam(ACW_50_mean ~ s(age, k = 3) + s(age, k=3, by=delayLength_order) + delayLength + sex, data = lobe_avg %>% filter(lobe == 'F'), random=list(lunaID=~1))
summary(gam.model)

draw(gam.model)

lunaize(
  ggplot(data = lobe_avg, 
         aes(x = age, y = ACW_50_mean, color=delayLength)) + geom_point() +
    geom_smooth(aes(group=delayLength), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
    labs(x = "", y = "ACW", fill = "") +
    theme_minimal()
) + 
  theme(
    legend.position = "left",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~lobe)




#avg the delay periods together

delay_avg <- allSub_delayPeriods_outlier %>%
  group_by(SubjectID, orig_parcelname, lobe, age, sex, visitno) %>%  # keep these in group
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    .groups = "drop"
  )%>% separate(SubjectID, c('lunaID','vdate'))

write.csv(delay_avg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam_variability.csv', row.names = F)


