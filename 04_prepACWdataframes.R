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

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# Age dataframes ----
merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$Subject <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")

# Delay 10 ----
acw_delay10 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay10.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_delay10_outlier <- acw_delay10 %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))

acw_delay10_outlier_avgTrials <- acw_delay10_outlier %>%
  group_by(lunaid, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

acw_delay10_outlier_avgTrials <- merge(chanLocs, acw_delay10_outlier_avgTrials, by ="channel")

acw_delay10_outlier_avgTrials <- acw_delay10_outlier_avgTrials %>%
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

acw_delay10_outlier_avgTrials <- acw_delay10_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_delay10_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_allChans_delay10.csv', row.names = F)

## Average across all regions and secondsEpochs ----
delay10_region_outlier_avgRegions <- acw_delay10_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(delay10_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_regionsAvgEpoch_delay10.csv', row.names = F)


## Average across all trials and channels ----
acw_delay10_outlier_avgChans <- acw_delay10_outlier %>%
  group_by(lunaid, secondsEpoch, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay10_outlier_avgChans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_avgChans_delay10.csv', row.names = F)

## Average across all trials and channels ----
acw_delay10_outlier_avgepoch <- acw_delay10_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))


## Average across all trials, channels, and secondsEpochs ----
acw_delay10_outlier_allAVG <- acw_delay10_outlier_avgChans %>%
  group_by(lunaid, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay10_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrialsChansEpochs_delay10.csv', row.names = F)

## Avg person ----
delay10_ageGroup_person <- acw_delay10_outlier_avgTrials %>%
  group_by(lunaid, secondsEpoch, region, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay10_ageGroup_person_outlier <- delay10_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(delay10_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgRegions_delay10.csv', row.names = F)



# Delay 8 ----
acw_delay8 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay8.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_delay8_outlier <- acw_delay8 %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))

acw_delay8_outlier_avgTrials <- acw_delay8_outlier %>%
  group_by(lunaid, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

acw_delay8_outlier_avgTrials <- merge(chanLocs, acw_delay8_outlier_avgTrials, by ="channel")

acw_delay8_outlier_avgTrials <- acw_delay8_outlier_avgTrials %>%
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

acw_delay8_outlier_avgTrials <- acw_delay8_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_delay8_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_allChans_delay8.csv', row.names = F)

## Average across all regions and secondsEpochs ----
delay8_region_outlier_avgRegions <- acw_delay8_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(delay8_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_regionsAvgEpoch_delay8.csv', row.names = F)

## Avg second epochs
acw_delay8_outlier_avgepoch <- acw_delay8_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

## Average across all trials and channels ----
acw_delay8_outlier_avgChans <- acw_delay8_outlier %>%
  group_by(lunaid, secondsEpoch, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay8_outlier_avgChans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_avgChans_delay8.csv', row.names = F)

## Average across all trials, channels, and secondsEpochs ----
acw_delay8_outlier_allAVG <- acw_delay8_outlier_avgChans %>%
  group_by(lunaid, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay8_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrialsChansEpochs_delay8.csv', row.names = F)

## Avg person ----
delay8_ageGroup_person <- acw_delay8_outlier_avgTrials %>%
  group_by(lunaid, secondsEpoch, region, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay8_ageGroup_person_outlier <- delay8_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(delay8_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgRegions_delay8.csv', row.names = F)


# Delay 6 ----
acw_delay6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay6.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_delay6_outlier <- acw_delay6 %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))

acw_delay6_outlier_avgTrials <- acw_delay6_outlier %>%
  group_by(lunaid, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

acw_delay6_outlier_avgTrials <- merge(chanLocs, acw_delay6_outlier_avgTrials, by ="channel")

acw_delay6_outlier_avgTrials <- acw_delay6_outlier_avgTrials %>%
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

acw_delay6_outlier_avgTrials <- acw_delay6_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_delay6_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_allChans_delay6.csv', row.names = F)

## Average across all regions and secondsEpochs ----
delay6_region_outlier_avgRegions <- acw_delay6_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(delay6_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_regionsAvgEpoch_delay6.csv', row.names = F)


## Average across all trials and channels ----
acw_delay6_outlier_avgChans <- acw_delay6_outlier %>%
  group_by(lunaid, secondsEpoch, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay6_outlier_avgChans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrials_avgChans_delay6.csv', row.names = F)

## Avg second epochs
acw_delay6_outlier_avgepoch <- acw_delay6_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

## Average across all trials, channels, and secondsEpochs ----
acw_delay6_outlier_allAVG <- acw_delay6_outlier_avgChans %>%
  group_by(lunaid, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_delay6_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgTrialsChansEpochs_delay6.csv', row.names = F)

## Avg person ----
delay6_ageGroup_person <- acw_delay6_outlier_avgTrials %>%
  group_by(lunaid, secondsEpoch, region, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay6_ageGroup_person_outlier <- delay6_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(delay6_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_avgRegions_delay6.csv', row.names = F)

# Combine Delay Periods ----
acw_delay10_outlier_avgepoch$epoch <- 'Delay10'
acw_delay8_outlier_avgepoch$epoch <- 'Delay8'
acw_delay6_outlier_avgepoch$epoch <- 'Delay6'

avgACW_allDelays <- rbind(acw_delay10_outlier_avgepoch, acw_delay8_outlier_avgepoch) %>% rbind(., acw_delay6_outlier_avgepoch)

### Average across delay epochs ----
avgACW_allDelaysavg <- avgACW_allDelays %>%
  group_by(lunaid, age, sex, visitno, channel) %>%
  summarize(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(avgACW_allDelaysavg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChannels_avgDelay.csv', row.names = F)



# Fix ----
acw_fix <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_allChans_acw_fix.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_fix_outlier <- acw_fix %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

acw_fix_outlier_avgTrials <- acw_fix_outlier %>%
  group_by(lunaid, Channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(acw_fix_outlier_avgTrials)[names(acw_fix_outlier_avgTrials) == "Channel"] <- "labels"

acw_fix_outlier_avgTrials <- merge(chanLocs, acw_fix_outlier_avgTrials, by ="labels")

acw_fix_outlier_avgTrials <- acw_fix_outlier_avgTrials %>%
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

acw_fix_outlier_avgTrials <- acw_fix_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_fix_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgTrials_allChans_fix.csv', row.names = F)

## Average across all regions and secondsEpochs ----
fix_region_outlier_avgRegions <- acw_fix_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(fix_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_regionsAvgEpoch_fix.csv', row.names = F)



## Average across all trials and channels ----
acw_fix_outlier_avgChans <- acw_fix_outlier %>%
  group_by(lunaid, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_fix_outlier_avgChans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgTrials_avgChans_fix.csv', row.names = F)

## Average across all trials, channels, and secondsEpochs ----
acw_fix_outlier_allAVG <- acw_fix_outlier_avgChans %>%
  group_by(lunaid, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_fix_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgTrialsChansEpochs_fix.csv', row.names = F)

## Avg person ----
fix_ageGroup_person <- acw_fix_outlier_avgTrials %>%
  group_by(lunaid, region, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

fix_ageGroup_person_outlier <- fix_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(fix_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgRegions_fix.csv', row.names = F)



# Rest ----
acw_eyesOpen <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesOpen.csv') %>%
  merge(., ageValues, by = "Subject")

acw_eyesOpen$epoch <- 'restEyesOpen'

acw_eyesClosed <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesClosed.csv') %>%
  merge(., ageValues, by = "Subject")

acw_eyesClosed$epoch <- 'restEyesClosed'

awc_rest <- rbind(acw_eyesOpen, acw_eyesClosed)

## Outlier detection ----
acw_rest_outlier <- awc_rest %>% group_by(Subject, epoch) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

acw_rest_outlier_avgTrials <- acw_rest_outlier %>%
  group_by(lunaid, Channel, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(acw_rest_outlier_avgTrials)[names(acw_rest_outlier_avgTrials) == "Channel"] <- "labels"

acw_rest_outlier_avgTrials <- merge(chanLocs, acw_rest_outlier_avgTrials, by ="labels")

acw_rest_outlier_avgTrials <- acw_rest_outlier_avgTrials %>%
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

acw_rest_outlier_avgTrials <- acw_rest_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_rest_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgTrials_allChans_rest.csv', row.names = F)

## Average across all regions and secondsEpochs ----
rest_region_outlier_avgRegions <- acw_rest_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(rest_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_regionsAvgEpoch_rest.csv', row.names = F)


## Average across all trials and channels ----
acw_rest_outlier_avgChans <- acw_rest_outlier %>%
  group_by(lunaid, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_rest_outlier_avgChans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgTrials_avgChans_rest.csv', row.names = F)

## Average across all trials, channels, and secondsEpochs ----
acw_rest_outlier_allAVG <- acw_rest_outlier_avgChans %>%
  group_by(lunaid, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_rest_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgTrialsChansEpochs_rest.csv', row.names = F)

## Avg person ----
rest_ageGroup_person <- acw_rest_outlier_avgTrials %>%
  group_by(lunaid, region, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

rest_ageGroup_person_outlier <- rest_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(rest_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgRegions_rest.csv', row.names = F)




acw_eyesOpen <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesOpen.csv') %>%
  merge(., ageValues, by = "Subject")

acw_eyesOpen$epoch <- 'restEyesOpen'

acw_eyesClosed <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesClosed.csv') %>%
  merge(., ageValues, by = "Subject")

acw_eyesClosed$epoch <- 'restEyesClosed'

awc_rest <- rbind(acw_eyesOpen, acw_eyesClosed)

## 2 second segments ----
acw <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/allSubjects_allChans_acw_2seconds_rest.csv') %>%
  merge(., ageValues, by = "Subject")

### Outlier detection ----
acw_rest_outlier <- acw %>% group_by(Subject, epoch) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

acw_rest_outlier_avgTrials <- acw_rest_outlier %>%
  group_by(lunaid, channel, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

## Average across all trials, channels, and secondsEpochs ----
acw_rest_outlier_allAVG <- acw_rest_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(acw_rest_outlier_allAVG, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/allSubjects_avgTrialsChansEpochs_rest.csv', row.names = F)


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

acw_rest_outlier_avgTrials <- merge(chanLocs, acw_rest_outlier_avgTrials, by ="channel")

acw_rest_outlier_avgTrials <- acw_rest_outlier_avgTrials %>%
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

acw_rest_outlier_avgTrials <- acw_rest_outlier_avgTrials %>% filter(region != "Temporal")

write.csv(acw_rest_outlier_avgTrials, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/allSubjects_avgTrials_allChans_rest.csv', row.names = F)

## Average across all regions and secondsEpochs ----
rest_region_outlier_avgRegions <- acw_rest_outlier_avgTrials %>%
  group_by(lunaid, age, sex, visitno, region, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% filter(region != 'Temporal')

write.csv(rest_region_outlier_avgRegions, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/allSubjects_regionsAvgEpoch_rest.csv', row.names = F)

## Avg person ----
rest_ageGroup_person <- acw_rest_outlier_avgTrials %>%
  group_by(lunaid, region, age, sex, visitno, epoch) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

rest_ageGroup_person_outlier <- rest_ageGroup_person %>% group_by(lunaid) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

write.csv(rest_ageGroup_person_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgRegions_rest.csv', row.names = F)



# Task state and Resting state all channels ----
avgACW_allDelaysavg <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChannels_avgDelay.csv')
avgACW_allDelaysavg$epoch <- 'delays'

acw_fix_outlier_avgTrials <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgTrials_allChans_fix.csv')
acw_fix_outlier_avgTrials$epoch <- 'fix'

acw_fix_outlier_avgTrials <- acw_fix_outlier_avgTrials %>%
  rename(channel = urchan)

delay_fix <- rbind(avgACW_allDelaysavg, acw_fix_outlier_avgTrials %>% select(c('lunaid', 'age', 'sex', 'visitno', 'channel', 'ACW_50', "ageGroup", 'epoch')))

task <- delay_fix %>%
  group_by(lunaid, age, sex, visitno, channel) %>%
  summarize(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% mutate(epoch='task')


acw_rest_outlier_avgTrials <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgTrials_allChans_rest.csv') %>% 
  filter(epoch =='restEyesOpen') %>%
  rename(channel = urchan) %>%
  mutate(epoch = 'rest') %>% select(c('lunaid', 'age', 'sex', 'visitno', 'channel', 'ACW_50', "ageGroup", 'epoch'))

rest_task <- rbind(task, acw_rest_outlier_avgTrials)


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv') %>% select(c("urchan", "labels"))
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

rest_task <- merge(chanLocs, rest_task, by ="channel")

rest_task <- rest_task %>%
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

write.csv(rest_task, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_taskRest_allChannels.csv', row.names = F)

# All Epochs all channels ----
avgACW_allDelaysavg <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChannels_avgDelay.csv')
avgACW_allDelaysavg$epoch <- 'delays'

acw_fix_outlier_avgTrials <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_avgTrials_allChans_fix.csv')
acw_fix_outlier_avgTrials$epoch <- 'fix'

acw_fix_outlier_avgTrials <- acw_fix_outlier_avgTrials %>%
  rename(channel = urchan)

delay_fix <- rbind(avgACW_allDelaysavg, acw_fix_outlier_avgTrials %>% select(c('lunaid', 'age', 'sex', 'visitno', 'channel', 'ACW_50', "ageGroup", 'epoch')))

task <- delay_fix %>%
  group_by(lunaid, age, sex, visitno, channel) %>%
  summarize(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% mutate(epoch='task')


acw_rest_outlier_avgTrials <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_avgTrials_allChans_rest.csv') %>% 
  rename(channel = urchan) %>% select(c('lunaid', 'age', 'sex', 'visitno', 'channel', 'ACW_50', "ageGroup", 'epoch'))

rest_task <- rbind(task, acw_rest_outlier_avgTrials)


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv') %>% select(c("urchan", "labels"))
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

rest_task <- merge(chanLocs, rest_task, by ="channel")

rest_task <- rest_task %>%
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

write.csv(rest_task, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/allSubjects_allEpochs_allChannels.csv', row.names = F)

