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

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","rest.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$SubjectID <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")
names(ageValues)[names(ageValues) == "visitDate"] <- "vdate"

# Task ----
## Delay ----
glasserINTs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam.csv')

glasserINTs <- glasserINTs %>%
  mutate(subjectID = paste(lunaID, vdate, sep = "_"))

glasserINTs <- glasserINTs %>% group_by(subjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTs <- glasserINTs %>% separate(subjectID, c('lunaid','vdate')) %>% select(-lunaID)
glasserINTs$epoch <- 'Delays'



## avg lobe
glasserINTs_avgLobe <-  glasserINTs %>% group_by(subjectID, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTs_avgLobe <- glasserINTs_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


glasserINTs_avgLobe <- glasserINTs_avgLobe %>% separate(subjectID, c('lunaid','vdate')) %>% filter(lobe != 'T')
glasserINTs_avgLobe$epoch <- 'Delays'

## Fix ----

fix_ints <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_fix.csv')

fix_ints <- fix_ints %>% group_by(SubjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier)) 

fix_ints <- fix_ints %>% separate(SubjectID, c('lunaid','vdate'))

fix_ints$epoch <- 'Fix'

fix_ints <- fix_ints %>% ungroup() %>% select(-c('ACW_0', 'OrigParcelName', "id", "label", "num", 'visitDate'))

## avg lobe
fix_ints_avgLobe <-  fix_ints %>% group_by(SubjectID, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

fix_ints_avgLobe <- fix_ints_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


fix_ints_avgLobe <- fix_ints_avgLobe %>% separate(SubjectID, c('lunaid','vdate')) %>% filter(lobe != 'T')

fix_ints_avgLobe$epoch <- 'Fix'



# Rest ----

glasserINTsRest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beam_rest_epochs.csv')

glasserINTsRest <- glasserINTsRest %>% group_by(SubjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTsRest <- glasserINTsRest %>% separate(SubjectID, c('lunaid','vdate'))

glasserINTsRest$epoch <- 'Rest'

glasserINTsRest <- glasserINTsRest %>% ungroup() %>% select(-c('ACW_0', 'OrigParcelName', "id", "label", "num", 'visitDate'))


## avg lobe 
glasserINTsRest_avgLobe <-  glasserINTsRest %>% group_by(lunaid, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTsRest_avgLobe <- glasserINTsRest_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


glasserINTsRest_avgLobe <- glasserINTsRest_avgLobe 

glasserINTsRest_avgLobe$epoch <- 'Rest'



# Combine all - all regions ---


allActivity <- rbind(glasserINTs, fix_ints) %>% rbind(., glasserINTsRest) 

write.csv(allActivity, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_allEpochs_allregions.csv', row.names = F)



# Combine all avg lobes ---


merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","rest.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$SubjectID <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")
names(ageValues)[names(ageValues) == "visitDate"] <- "vdate"

allActivity <- rbind(glasserINTs_avgLobe, fix_ints_avgLobe) %>% rbind(., glasserINTsRest_avgLobe) %>% ungroup() %>% select(-lobe) %>% merge(., ageValues, by = c('lunaid', 'vdate', 'visitno'))

write.csv(allActivity, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_allEpochs.csv', row.names = F)



# Variability ----
# Task ----
## Delay ----
glasserINTs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_MGS_avgDelayLengths_beam_variability.csv')

glasserINTs <- glasserINTs %>%
  mutate(subjectID = paste(lunaID, vdate, sep = "_"))

glasserINTs <- glasserINTs %>% group_by(subjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTs$epoch <- 'Delays'



## avg lobe
glasserINTs_avgLobe <-  glasserINTs %>% group_by(subjectID, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTs_avgLobe <- glasserINTs_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


glasserINTs_avgLobe <- glasserINTs_avgLobe %>% separate(subjectID, c('lunaid','vdate')) %>% filter(lobe != 'T')
glasserINTs_avgLobe$epoch <- 'Delays'

## Fix ----

fix_ints <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_fix_variability.csv')

fix_ints <- fix_ints %>% group_by(SubjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier)) 


fix_ints$epoch <- 'Fix'

fix_ints <- fix_ints %>% ungroup() %>% select(-c('ACW_0', 'OrigParcelName', "id", "label", "num", 'visitDate'))

## avg lobe
fix_ints_avgLobe <-  fix_ints %>% group_by(SubjectID, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

fix_ints_avgLobe <- fix_ints_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


fix_ints_avgLobe <- fix_ints_avgLobe %>% separate(SubjectID, c('lunaid','vdate')) %>% filter(lobe != 'T')

fix_ints_avgLobe$epoch <- 'Fix'



# Rest ----

glasserINTsRest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beam_rest_epochs_variability.csv')

glasserINTsRest <- glasserINTsRest %>% group_by(SubjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))


glasserINTsRest$epoch <- 'Rest'

glasserINTsRest <- glasserINTsRest %>% ungroup() %>% select(-c('ACW_0', 'OrigParcelName', "id", "label", "num", 'visitDate'))


## avg lobe 
glasserINTsRest_avgLobe <-  glasserINTsRest %>% group_by(lunaid, lobe, visitno) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTsRest_avgLobe <- glasserINTsRest_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal",
                         "T" = "Temporal"))


glasserINTsRest_avgLobe <- glasserINTsRest_avgLobe 

glasserINTsRest_avgLobe$epoch <- 'Rest'



# Combine all - all regions ---
fix_ints <- fix_ints %>% separate(SubjectID, c('lunaid','vdate'))
glasserINTsRest <- glasserINTsRest %>% separate(SubjectID, c('lunaid','vdate'))
glasserINTs <- glasserINTs %>% separate(subjectID, c('lunaid','vdate')) %>% select(-lunaID)

allActivity <- rbind(glasserINTs%>% select(-lunaID),  fix_ints) %>% rbind(., glasserINTsRest) 

write.csv(allActivity, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_allEpochs_allregions_variability.csv', row.names = F)



# Combine all avg lobes ---


merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","rest.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$SubjectID <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")
names(ageValues)[names(ageValues) == "visitDate"] <- "vdate"

allActivity <- rbind(glasserINTs_avgLobe, fix_ints_avgLobe) %>% rbind(., glasserINTsRest_avgLobe) %>% ungroup() %>% select(-lobe) %>% merge(., ageValues, by = c('lunaid', 'vdate', 'visitno'))

write.csv(allActivity, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_allEpochs.csv', row.names = F)
