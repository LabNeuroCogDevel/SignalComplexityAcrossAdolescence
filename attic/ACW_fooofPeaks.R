
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

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()

# Age dataframes ----
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age", "sex")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/taskRest.csv') %>% select(-X)

names(task_rest)[names(task_rest) == "lunaid"] <- "lunaID"


# Load in peak measures ----
peaks <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/allSubjectsPeakMeasures_20230710.csv') %>% filter(Condition == 'eyesOpen') %>% select(-X)

peaks <- peaks %>%
  mutate(Band = case_when(
    Center.Frequency >= 0.5  & Center.Frequency < 4   ~ "Delta",
    Center.Frequency >= 4    & Center.Frequency < 8   ~ "Theta",
    Center.Frequency >= 8    & Center.Frequency < 12  ~ "Alpha",
    Center.Frequency >= 12   & Center.Frequency < 30  ~ "Beta",
    Center.Frequency >= 30   & Center.Frequency < 100 ~ "Gamma",
    TRUE ~ "Other"
  )) %>% filter(Band != "Delta")

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/resources/ChannelLocs.csv')%>% select(-type)
names(chanLocs)[names(chanLocs) == "urchan"] <- "Channel"

peaks_chanLocs <- merge(chanLocs, peaks, by = 'Channel')

peaks_chanLocs <- peaks_chanLocs %>%
  mutate(region = case_when(
    startsWith(labels, "F") ~ "Frontal",
    startsWith(labels, "A") ~ "Frontal",
    startsWith(labels, "Cz") | startsWith(labels, "C1")  | startsWith(labels, "C2") ~ "Frontal",
    startsWith(labels, "C3") | startsWith(labels, "C4")  | startsWith(labels, "C5") | startsWith(labels, "C6") | startsWith(labels, "CPz") ~ "SensoryMotor",
    startsWith(labels, "CP") ~ "SensoryMotor",
    startsWith(labels, "P") ~ "Parietal",
    startsWith(labels, "O") ~ "Occipital",
    startsWith(labels, "I") ~ "Occipital",
    startsWith(labels, "T") ~ "Temporal"
  ))

peaks_chanLocs <- peaks_chanLocs %>% filter(region != "Temporal") %>% filter(region != "SensoryMotor")

peaks_chanLocs_outlier <- peaks_chanLocs %>% group_by(Subject, region, Band) %>%
  mutate(across(c('Power'), naoutlier)) 

peaks_chanLocs_outlier_avgRegions <- peaks_chanLocs_outlier %>%
  group_by(Subject, region, Band) %>%
  summarize(
    Power = mean(Power, na.rm = TRUE))

peaks_chanLocs_outlier_avgRegions_outlier <- peaks_chanLocs_outlier_avgRegions %>% group_by(Band) %>%
  mutate(across(c('Power'), naoutlier)) 


peaks_chanLocs_outlier_avgRegions_outlier <- merge(peaks_chanLocs_outlier_avgRegions_outlier, ageValues, by = c('Subject')) %>% separate(Subject, c('lunaID','vdate'))

# Merge peaks with ACW ----

ACW_peaks <- merge(peaks_chanLocs_outlier_avgRegions_outlier, task_rest %>% filter(epoch == 'restEyesOpen'), by = c('lunaID', 'visitno', 'age', 'region', 'sex'))

## Theta ----
ACW_peaks_Theta <- ACW_peaks %>% filter(Band == 'Theta')


lunaize(ggplot(data = ACW_peaks_Theta, aes(y = ACW_50, x = Power, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='lm', alpha=0.5, linewidth = 1.5)+ xlab('theta power') +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)


### Frontal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Theta %>% filter(region=='Frontal'))
summary(lm.model)

### Parietal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Theta %>% filter(region=='Parietal'))
summary(lm.model)

### Occipital ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Theta %>% filter(region=='Occipital'))
summary(lm.model)


### Age Interaction ----

lm.interactionmodel <- lmerTest::lmer(ACW_50 ~ Power*age + region + (1|lunaID), data = ACW_peaks_Theta)
summary(lm.interactionmodel)

ggpredictinteractiondata <- (ggpredict(lm.interactionmodel, terms = c("Power", "age", 'region')))

lunaize(ggplot( ) + 
          geom_point(data = ACW_peaks_Theta, aes(y = ACW_50, x = Power, shape=region, color = region), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, color=group), data=ggpredictinteractiondata, size = 2)+ 
          theme_minimal()+ 
          labs(y = "ACW", x = "theta Power")) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2))+ facet_wrap(~region)



## Alpha ----
ACW_peaks_Alpha <- ACW_peaks %>% filter(Band == 'Alpha')


lunaize(ggplot(data = ACW_peaks_Alpha, aes(y = Power, x = age, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='gam', alpha=0.5, linewidth = 1.5) +xlab('')+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)

lunaize(ggplot(data = ACW_peaks_Alpha, aes(y = ACW_50, x = Power, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='lm', alpha=0.5, linewidth = 1.5) +xlab('alpha power')+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)


### Frontal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Alpha %>% filter(region=='Frontal'))
summary(lm.model)

ACW_peaks_Alpha %>% 
  filter(region == 'Frontal') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))

### Parietal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Alpha %>% filter(region=='Parietal'))
summary(lm.model)
p.adjust(0.007, method = 'bonferroni', n=3)

ACW_peaks_Alpha %>% 
  filter(region == 'Parietal') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))


### Occipital ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Alpha %>% filter(region=='Occipital'))
summary(lm.model)

ACW_peaks_Alpha %>% 
  filter(region == 'Occipital') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))

### Age Interaction ----

lm.interactionmodel <- lmerTest::lmer(ACW_50 ~ Power*age + region + (1|lunaID), data = ACW_peaks_Alpha)
summary(lm.interactionmodel)
  
ggpredictinteractiondata <- (ggpredict(lm.interactionmodel, terms = c("Power", "age", 'region')))

lunaize(ggplot( ) + 
          geom_point(data = ACW_peaks_Alpha, aes(y = ACW_50, x = Power, shape=region, color = region), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, color=group), data=ggpredictinteractiondata, size = 2)+ 
          theme_minimal()+ 
          labs(y = "ACW", x = "Alpha Power")) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2)) + facet_wrap(~region)


## Beta ----
ACW_peaks_Beta <- ACW_peaks %>% filter(Band == 'Beta')

lunaize(ggplot(data = ACW_peaks_Beta, aes(y = Power, x = age, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='lm', alpha=0.5, linewidth = 1.5) +xlab('')+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)

lunaize(ggplot(data = ACW_peaks_Beta, aes(y = ACW_50, x = Power, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='lm', alpha=0.5, linewidth = 1.5)+ xlab('beta power') +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)


### Frontal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Beta %>% filter(region=='Frontal'))
summary(lm.model)

ACW_peaks_Beta %>% 
  filter(region == 'Frontal') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))


### Parietal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Beta %>% filter(region=='Parietal'))
summary(lm.model)

ACW_peaks_Beta %>% 
  filter(region == 'Parietal') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))


### Occipital ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Beta %>% filter(region=='Occipital'))
summary(lm.model)

ACW_peaks_Beta %>% 
  filter(region == 'Occipital') %>% 
  summarise(correlation = cor(ACW_50, Power, use = 'complete.obs'))

### Age Interaction ----

lm.interactionmodel <- lmerTest::lmer(ACW_50 ~ Power*age + region + (1|lunaID), data = ACW_peaks_Beta)
summary(lm.interactionmodel)

ggpredictinteractiondata <- (ggpredict(lm.interactionmodel, terms = c("Power", "age", 'region')))

lunaize(ggplot( ) + 
          geom_point(data = ACW_peaks_Beta, aes(y = ACW_50, x = Power, shape=region, color = region), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, color=group), data=ggpredictinteractiondata, size = 2)+ 
          theme_minimal()+ 
          labs(y = "ACW", x = "Power")) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2)) + facet_wrap(~region)



## Gamma ----
ACW_peaks_Gamma <- ACW_peaks %>% filter(Band == 'Gamma')

lunaize(ggplot(data = ACW_peaks_Gamma, aes(y = ACW_50, x = Power, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaID, Band)), alpha = 0.2) +
  geom_smooth(aes(group=c(Band)), method='lm', alpha=0.5, linewidth = 1.5)+ xlab('gamma power')+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)


### Frontal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Gamma %>% filter(region=='Frontal'))
summary(lm.model)

### Parietal ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Gamma %>% filter(region=='Parietal'))
summary(lm.model)

### Occipital ----
lm.model <- lmerTest::lmer(ACW_50 ~ Power + age + (1|lunaID), data = ACW_peaks_Gamma %>% filter(region=='Occipital'))
summary(lm.model)


### Age Interaction ----

lm.interactionmodel <- lmerTest::lmer(ACW_50 ~ Power*age + region + (1|lunaID), data = ACW_peaks_Gamma)
summary(lm.interactionmodel)

ggpredictinteractiondata <- (ggpredict(lm.interactionmodel, terms = c("Power", "age", 'region')))

lunaize(ggplot( ) + 
          geom_point(data = ACW_peaks_Gamma, aes(y = ACW_50, x = Power, shape=region, color = region), alpha = 0.3) + 
          geom_smooth(aes(x = x, y = predicted, group=group, color=group), data=ggpredictinteractiondata, size = 2)+ 
          theme_minimal()+ 
          labs(y = "ACW", x = "Power")) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20,  color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2)) + facet_wrap(~region)


