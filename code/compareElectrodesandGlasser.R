
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


avgACW_allEpochs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allEpochs_regions.csv') %>% 
  filter(epoch == 'Delays') %>%
  select(-ACW_0)

avgACW_allEpochs <- avgACW_allEpochs %>% group_by(region) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam.csv')
glasserINTs <- glasserINTs %>%
  mutate(subjectID = paste(lunaID, vdate, sep = "_"))

glasserINTs <- glasserINTs %>% group_by(lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))


glasserINTs_avgLobe <-  glasserINTs %>% group_by(subjectID, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTs_avgLobe <- glasserINTs_avgLobe %>%
  mutate(region = recode(lobe,
                       "F" = "Frontal",
                       "O" = "Occipital",
                       "P" = "Parietal"))


glasserINTs_avgLobe <- glasserINTs_avgLobe %>% separate(subjectID, c('lunaid','vdate'))

glasser_electrode <- merge(glasserINTs_avgLobe, avgACW_allEpochs, by = c('lunaid', 'visitno', 'region'))
glasser_electrode <- glasser_electrode %>%
  rename(INT_electrode = ACW_50) %>% select(-lobe)


(lunaize(ggplot(data = glasser_electrode, aes(x = INT_glasser, y = INT_electrode, color=region))) + geom_point(alpha =0.3) +
          geom_smooth(aes(group=1), method="lm", alpha=0.3, linewidth = 1.5)+ facet_wrap(~region)+
          scale_color_manual(values = c("#012b0d", "#90b8a6", "#27763d")) + 
          theme(
            legend.position = "none",
            axis.text = element_text(size = 20, color = c("black")),
            axis.title.x = element_text(size = 20, color = c("black")),
            axis.title.y = element_text(size = 20, color = c("black")),
            axis.line = element_line(linewidth = .2), 
            axis.ticks = element_line(linewidth = .2),
            strip.background = element_rect('white')))

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Supp/electrodeVsGlasser.pdf", device = "pdf", dpi = 500, width = 10, height = 5)


lm.model <- lme4::lmer(INT_glasser ~ INT_electrode + region + (1|lunaid), data = glasser_electrode)
summary(lm.model)

glasser_electrode_frontal <- glasser_electrode %>% filter(region == 'Frontal')
cor(glasser_electrode_frontal$INT_glasser, glasser_electrode_frontal$INT_electrode, use = "pairwise.complete.obs")

lm.model <- lmerTest::lmer(INT_glasser ~ INT_electrode + (1|lunaid), data = glasser_electrode_frontal)
summary(lm.model)


glasser_electrode_parietal <- glasser_electrode %>% filter(region == 'Parietal')
cor(glasser_electrode_parietal$INT_glasser, glasser_electrode_parietal$INT_electrode, use = "pairwise.complete.obs")

lm.model <- lmerTest::lmer(INT_glasser ~ INT_electrode + (1|lunaid), data = glasser_electrode_parietal)
summary(lm.model)

glasser_electrode_occipital <- glasser_electrode %>% filter(region == 'Occipital')
cor(glasser_electrode_occipital$INT_glasser, glasser_electrode_occipital$INT_electrode,  use = "pairwise.complete.obs")

lm.model <- lmerTest::lmer(INT_glasser ~ INT_electrode + (1|lunaid), data = glasser_electrode_occipital)
summary(lm.model)


## channel level 

avgACW_allDelaysavg <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChannels_avgDelay.csv')
dlpfc <- avgACW_allDelaysavg %>% filter(channel == 6)

dlpfc <- dlpfc %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam.csv')
glasserINTs <- glasserINTs %>%
  mutate(subjectID = paste(lunaID, vdate, sep = "_"))
glasserINTs_dlpfc <- glasserINTs %>% filter(orig_parcelname == 'L_9-46d_ROI')

glasserINTs_dlpfc <- glasserINTs_dlpfc %>%
  rename(lunaid = lunaID)

glasserINTs_dlpfc <- glasserINTs_dlpfc %>% group_by(lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasser_electrode_dlpfc <- merge(glasserINTs_dlpfc, dlpfc, by = c('lunaid', 'visitno', 'age', 'sex'), suffixes = c('_glasser', '_electrode'))

(lunaize(ggplot(data = glasser_electrode_dlpfc, aes(x = ACW_50_glasser, y = ACW_50_electrode))) + geom_point(alpha =0.3) +
    geom_smooth(aes(group=1), method="lm", alpha=0.3, linewidth = 1.5) +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, color = c("black")),
      axis.title.x = element_text(size = 20, color = c("black")),
      axis.title.y = element_text(size = 20, color = c("black")),
      axis.line = element_line(linewidth = .2), 
      axis.ticks = element_line(linewidth = .2),
      strip.background = element_rect('white')))

cor(glasser_electrode_dlpfc$ACW_50_glasser,
    glasser_electrode_dlpfc$ACW_50_electrode,
    use = "complete.obs")
