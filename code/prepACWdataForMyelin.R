
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
library(eegkit)


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


acw_delay10_outlier_avgTrials_avgEpochs <- acw_delay10_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

acw_delay10_outlier_avgTrials_avgEpochs$epoch <- 'delay10'

# Delay 8 ----
acw_delay8 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay8.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_delay8_outlier <- acw_delay8 %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


acw_delay8_outlier_avgTrials_avgEpochs <- acw_delay8_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))


acw_delay8_outlier_avgTrials_avgEpochs$epoch <- 'delay8'


# Delay 6 ----
acw_delay6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay6.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_delay6_outlier <- acw_delay6 %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


acw_delay6_outlier_avgTrials_avgEpochs <- acw_delay6_outlier %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))


acw_delay6_outlier_avgTrials_avgEpochs$epoch <- 'delay6'

# Combine Delay Epochs ----
acw_allDelays <- rbind(acw_delay10_outlier_avgTrials_avgEpochs, acw_delay8_outlier_avgTrials_avgEpochs) %>% rbind(., acw_delay6_outlier_avgTrials_avgEpochs)

acw_avgDelays <- acw_allDelays %>%
  group_by(lunaid, channel, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

acw_avgDelays_chanLocs <- merge(chanLocs %>% select("channel", "labels"), acw_avgDelays, by ="channel")



## Average acw measures across electrodes ----

acw_avgDelays_wide <- acw_avgDelays_chanLocs %>% select(-channel) %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(ACW_0, ACW_50),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like ACW_0_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

acw_avgDelays_wide <- combine_electrodes(df = acw_avgDelays_wide, electrode_list = c("F7_ACW_50", "F8_ACW_50"), output_name = "vlpfc_ACW_50")
acw_avgDelays_wide <- combine_electrodes(df = acw_avgDelays_wide, electrode_list = c("AF5_ACW_50", "AF6_ACW_50", "F3_ACW_50", "F4_ACW_50", "F5_ACW_50", "F6_ACW_50"), output_name = "dlpfc_ACW_50")
acw_avgDelays_wide <- combine_electrodes(df = acw_avgDelays_wide, electrode_list = c("AF1_ACW_50", "AF2_ACW_50", "F1_ACW_50", "F2_ACW_50"), output_name = "spfc_ACW_50")
acw_avgDelays_wide <- combine_electrodes(df = acw_avgDelays_wide, electrode_list = c("FC1_ACW_50", "FC2_ACW_50", "FC3_ACW_50", "FC4_ACW_50"), output_name = "motor_ACW_50")

saveRDS(acw_avgDelays_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/ACW_readyforMyelin.RDS")


## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

myelin.acw.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, acw_avgDelays_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.acw.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/R1_ACW_allDelays.RDS")

# Rest Eyes Open----
acw_eyesOpen <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesOpen.csv') %>%
  merge(., ageValues, by = "Subject")


## Outlier detection ----
acw_eyesOpen_outlier <- acw_eyesOpen %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(acw_eyesOpen_outlier)[names(acw_eyesOpen_outlier) == "Channel"] <- "labels"

acw_eyesOpen_outlier <- merge(chanLocs %>% select("labels"), acw_eyesOpen_outlier, by ="labels")

## Average acw measures across electrodes ----

acw_rest_outlier_wide <- acw_eyesOpen_outlier %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(ACW_0, ACW_50),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like ACW_0_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

acw_rest_outlier_wide <- combine_electrodes(df = acw_rest_outlier_wide, electrode_list = c("F7_ACW_50", "F8_ACW_50"), output_name = "vlpfc_ACW_50")
acw_rest_outlier_wide <- combine_electrodes(df = acw_rest_outlier_wide, electrode_list = c("AF5_ACW_50", "AF6_ACW_50", "F3_ACW_50", "F4_ACW_50", "F5_ACW_50", "F6_ACW_50"), output_name = "dlpfc_ACW_50")
acw_rest_outlier_wide <- combine_electrodes(df = acw_rest_outlier_wide, electrode_list = c("AF1_ACW_50", "AF2_ACW_50", "F1_ACW_50", "F2_ACW_50"), output_name = "spfc_ACW_50")
acw_rest_outlier_wide <- combine_electrodes(df = acw_rest_outlier_wide, electrode_list = c("FC1_ACW_50", "FC2_ACW_50", "FC3_ACW_50", "FC4_ACW_50"), output_name = "motor_ACW_50")

saveRDS(acw_rest_outlier_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/ACW_readyforMyelin_resteyesOpen.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

acw_rest_outlier_wide$lunaid <- as.integer(acw_rest_outlier_wide$lunaid)

myelin.acw.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, acw_rest_outlier_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.acw.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/R1_ACW_resteyesOpen.RDS")




# Rest Eyes Closed----
acw_eyesClosed <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesClosed.csv') %>%
  merge(., ageValues, by = "Subject")


## Outlier detection ----
acw_eyesClosed_outlier <- acw_eyesClosed %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(acw_eyesClosed_outlier)[names(acw_eyesClosed_outlier) == "Channel"] <- "labels"

acw_eyesClosed_outlier <- merge(chanLocs %>% select("labels"), acw_eyesClosed_outlier, by ="labels")

## Average acw measures across electrodes ----

acw_eyesClosed_outlier_wide <- acw_eyesClosed_outlier %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(ACW_0, ACW_50),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like ACW_0_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

acw_eyesClosed_outlier_wide <- combine_electrodes(df = acw_eyesClosed_outlier_wide, electrode_list = c("F7_ACW_50", "F8_ACW_50"), output_name = "vlpfc_ACW_50")
acw_eyesClosed_outlier_wide <- combine_electrodes(df = acw_eyesClosed_outlier_wide, electrode_list = c("AF5_ACW_50", "AF6_ACW_50", "F3_ACW_50", "F4_ACW_50", "F5_ACW_50", "F6_ACW_50"), output_name = "dlpfc_ACW_50")
acw_eyesClosed_outlier_wide <- combine_electrodes(df = acw_eyesClosed_outlier_wide, electrode_list = c("AF1_ACW_50", "AF2_ACW_50", "F1_ACW_50", "F2_ACW_50"), output_name = "spfc_ACW_50")
acw_eyesClosed_outlier_wide <- combine_electrodes(df = acw_eyesClosed_outlier_wide, electrode_list = c("FC1_ACW_50", "FC2_ACW_50", "FC3_ACW_50", "FC4_ACW_50"), output_name = "motor_ACW_50")

saveRDS(acw_eyesClosed_outlier_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/ACW_readyforMyelin_resteyesClosed.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

acw_eyesClosed_outlier_wide$lunaid <- as.integer(acw_eyesClosed_outlier_wide$lunaid)

myelin.acw.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, acw_eyesClosed_outlier_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.acw.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/R1_ACW_resteyesClosed.RDS")



# Fix ----
acw_fix <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_allChans_acw_fix.csv') %>%
  merge(., ageValues, by = "Subject")

## Outlier detection ----
acw_fix_outlier <- acw_fix %>% group_by(Subject) %>%
  mutate(across(c("ACW_0","ACW_50"), naoutlier)) %>% separate(Subject, c('lunaid','vdate'))

acw_fix_outlier_avgTrials <- acw_fix_outlier %>%
  group_by(lunaid, labels, age, sex, visitno) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(acw_fix_outlier)[names(acw_fix_outlier) == "Channel"] <- "labels"

acw_fix_outlier_avgTrials <- merge(chanLocs %>% select("labels"), acw_fix_outlier_avgTrials, by ="labels")

## Average acw measures across electrodes ----

acw_fix_outlier_avgTrials_wide <- acw_fix_outlier_avgTrials %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(ACW_0, ACW_50),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like ACW_0_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

acw_fix_outlier_avgTrials_wide <- combine_electrodes(df = acw_fix_outlier_avgTrials_wide, electrode_list = c("F7_ACW_50", "F8_ACW_50"), output_name = "vlpfc_ACW_50")
acw_fix_outlier_avgTrials_wide <- combine_electrodes(df = acw_fix_outlier_avgTrials_wide, electrode_list = c("AF5_ACW_50", "AF6_ACW_50", "F3_ACW_50", "F4_ACW_50", "F5_ACW_50", "F6_ACW_50"), output_name = "dlpfc_ACW_50")
acw_fix_outlier_avgTrials_wide <- combine_electrodes(df = acw_fix_outlier_avgTrials_wide, electrode_list = c("AF1_ACW_50", "AF2_ACW_50", "F1_ACW_50", "F2_ACW_50"), output_name = "spfc_ACW_50")
acw_fix_outlier_avgTrials_wide <- combine_electrodes(df = acw_fix_outlier_avgTrials_wide, electrode_list = c("FC1_ACW_50", "FC2_ACW_50", "FC3_ACW_50", "FC4_ACW_50"), output_name = "motor_ACW_50")

saveRDS(acw_fix_outlier_avgTrials_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/ACW_readyforMyelin_fix.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

acw_fix_outlier_avgTrials_wide$lunaid <- as.integer(acw_fix_outlier_avgTrials_wide$lunaid)

myelin.acw.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, acw_fix_outlier_avgTrials_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.acw.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/R1_ACW_fix.RDS")





