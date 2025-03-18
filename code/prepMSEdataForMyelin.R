
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
mse_delay10 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgChans_MSE20_delay10.csv')

delay10_avgTrials_chans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay10.csv')

  
delay10_secondsEpoch <- delay10_avgTrials_chans %>%
  group_by(lunaid, labels, age, sex, visitno) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay10_secondsEpoch_outlier <- delay10_secondsEpoch %>% group_by(lunaid) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

delay10_secondsEpoch_outlier$epoch <- 'delay10'

# Delay 8 ----
mse_delay8 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgChans_MSE20_delay8.csv')

delay8_avgTrials_chans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay8.csv')


delay8_secondsEpoch <- delay8_avgTrials_chans %>%
  group_by(lunaid, labels, age, sex, visitno) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay8_secondsEpoch_outlier <- delay8_secondsEpoch %>% group_by(lunaid) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

delay8_secondsEpoch_outlier$epoch <- 'delay8'

# Delay 6 ----
mse_delay6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgChans_MSE20_delay6.csv')

delay6_avgTrials_chans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay6.csv')


delay6_secondsEpoch <- delay6_avgTrials_chans %>%
  group_by(lunaid, labels, age, sex, visitno) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

delay6_secondsEpoch_outlier <- delay6_secondsEpoch %>% group_by(lunaid) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

delay6_secondsEpoch_outlier$epoch <- 'delay6'


# Combine Delay Epochs ----
mse_allDelays <- rbind(delay10_secondsEpoch_outlier, delay8_secondsEpoch_outlier) %>% rbind(., delay6_secondsEpoch_outlier)

mse_avgDelays <- mse_allDelays %>%
  group_by(lunaid, labels, age, sex, visitno) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

mse_avgDelays_outlier <- mse_avgDelays %>% group_by(lunaid) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

mse_avgDelays_chanLocs <- merge(chanLocs %>% select("channel", "labels"), mse_avgDelays_outlier, by ="labels")



## Average mse measures across electrodes ----

mse_avgDelays_wide <- mse_avgDelays_chanLocs %>% select(-channel) %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(Var1),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like Var1_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

mse_avgDelays_wide <- combine_electrodes(df = mse_avgDelays_wide, electrode_list = c("F7_Var1", "F8_Var1"), output_name = "vlpfc_Var1")
mse_avgDelays_wide <- combine_electrodes(df = mse_avgDelays_wide, electrode_list = c("AF5_Var1", "AF6_Var1", "F3_Var1", "F4_Var1", "F5_Var1", "F6_Var1"), output_name = "dlpfc_Var1")
mse_avgDelays_wide <- combine_electrodes(df = mse_avgDelays_wide, electrode_list = c("AF1_Var1", "AF2_Var1", "F1_Var1", "F2_Var1"), output_name = "spfc_Var1")
mse_avgDelays_wide <- combine_electrodes(df = mse_avgDelays_wide, electrode_list = c("FC1_Var1", "FC2_Var1", "FC3_Var1", "FC4_Var1"), output_name = "motor_Var1")

mse_avgDelays_wide <- mse_avgDelays_wide %>%
  mutate(across(c("vlpfc_Var1", "dlpfc_Var1", "spfc_Var1", "motor_Var1"), naoutlier)) %>% ungroup()

saveRDS(mse_avgDelays_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/mse_readyforMyelin.RDS")


## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

myelin.mse.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, mse_avgDelays_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.mse.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/R1_mse_allDelays.RDS")

# Rest Eyes Open----

rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/allSubjects_allChans_MSE20_2seconds_rest.csv') %>%
  merge(., ageValues, by="Subject") %>% separate(Subject, c('lunaid','vdate'))

rest_wholeBrain <- rest %>%
  group_by(lunaid, epoch, age, sex, visitno, channel) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

rest_wholeBrain_outlier <- rest_wholeBrain %>% group_by(lunaid,epoch) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

rest_wholeBrain_outlier$epoch <- as.factor(rest_wholeBrain_outlier$epoch)

mse_eyesOpen_outlier <- merge(chanLocs %>% select("channel", "labels"), rest_wholeBrain_outlier, by ="channel") %>% dplyr::filter(epoch=='restEyesOpen')


## Average mse measures across electrodes ----
mse_eyesOpen_outlier <- mse_eyesOpen_outlier %>% select(-channel)

mse_rest_outlier_wide <- mse_eyesOpen_outlier %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(Var1),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like Var1_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

mse_rest_outlier_wide <- combine_electrodes(df = mse_rest_outlier_wide, electrode_list = c("F7_Var1", "F8_Var1"), output_name = "vlpfc_Var1")
mse_rest_outlier_wide <- combine_electrodes(df = mse_rest_outlier_wide, electrode_list = c("AF5_Var1", "AF6_Var1", "F3_Var1", "F4_Var1", "F5_Var1", "F6_Var1"), output_name = "dlpfc_Var1")
mse_rest_outlier_wide <- combine_electrodes(df = mse_rest_outlier_wide, electrode_list = c("AF1_Var1", "AF2_Var1", "F1_Var1", "F2_Var1"), output_name = "spfc_Var1")
mse_rest_outlier_wide <- combine_electrodes(df = mse_rest_outlier_wide, electrode_list = c("FC1_Var1", "FC2_Var1", "FC3_Var1", "FC4_Var1"), output_name = "motor_Var1")

mse_rest_outlier_wide <- mse_rest_outlier_wide %>%
  mutate(across(c("vlpfc_Var1", "dlpfc_Var1", "spfc_Var1", "motor_Var1"), naoutlier)) %>% ungroup()


saveRDS(mse_rest_outlier_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/mse_readyforMyelin_resteyesOpen.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

mse_rest_outlier_wide$lunaid <- as.integer(mse_rest_outlier_wide$lunaid)

myelin.mse.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, mse_rest_outlier_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.mse.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/R1_mse_resteyesOpen.RDS")




# Rest Eyes Closed----
mse_eyesClosed_outlier <- merge(chanLocs %>% select("channel", "labels"), rest_wholeBrain_outlier, by ="channel") %>% dplyr::filter(epoch=='restEyesClosed')

## Average mse measures across electrodes ----
mse_eyesClosed_outlier <- mse_eyesClosed_outlier %>% select(-channel)

mse_eyesClosed_outlier_wide <- mse_eyesClosed_outlier %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(Var1),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like Var1_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

mse_eyesClosed_outlier_wide <- combine_electrodes(df = mse_eyesClosed_outlier_wide, electrode_list = c("F7_Var1", "F8_Var1"), output_name = "vlpfc_Var1")
mse_eyesClosed_outlier_wide <- combine_electrodes(df = mse_eyesClosed_outlier_wide, electrode_list = c("AF5_Var1", "AF6_Var1", "F3_Var1", "F4_Var1", "F5_Var1", "F6_Var1"), output_name = "dlpfc_Var1")
mse_eyesClosed_outlier_wide <- combine_electrodes(df = mse_eyesClosed_outlier_wide, electrode_list = c("AF1_Var1", "AF2_Var1", "F1_Var1", "F2_Var1"), output_name = "spfc_Var1")
mse_eyesClosed_outlier_wide <- combine_electrodes(df = mse_eyesClosed_outlier_wide, electrode_list = c("FC1_Var1", "FC2_Var1", "FC3_Var1", "FC4_Var1"), output_name = "motor_Var1")

mse_eyesClosed_outlier_wide <- mse_eyesClosed_outlier_wide %>%
  mutate(across(c("vlpfc_Var1", "dlpfc_Var1", "spfc_Var1", "motor_Var1"), naoutlier)) %>% ungroup()


saveRDS(mse_eyesClosed_outlier_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/mse_readyforMyelin_resteyesClosed.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

mse_eyesClosed_outlier_wide$lunaid <- as.integer(mse_eyesClosed_outlier_wide$lunaid)

myelin.mse.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, mse_eyesClosed_outlier_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.mse.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/R1_mse_resteyesClosed.RDS")



# Fix ----
fix <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/allSubjects_allChans_mse_fix.csv') %>%
  merge(., ageValues, by="Subject") %>% separate(Subject, c('lunaid','vdate'))

fix_chans <- fix %>%
  group_by(lunaid, age, sex, visitno, channel) %>%
  summarize(
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

fix_chans_outlier <- fix_chans %>% group_by(lunaid) %>%
  mutate(across(c("Var1"), naoutlier)) %>% ungroup()

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')
names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

fix_chans_outlier <- merge(chanLocs %>% select("labels", "channel"), fix_chans_outlier, by ="channel")

## Average mse measures across electrodes ----
fix_chans_outlier <- fix_chans_outlier %>% select(-channel)

mse_fix_outlier_avgTrials_wide <- fix_chans_outlier %>%
  pivot_wider(
    names_from = labels,  # Labels become column suffixes
    values_from = c(Var1),  # Values to be spread out
    names_glue = "{labels}_{.value}"  # Custom column names like Var1_label
  )

combine_electrodes <- function(df, electrode_list, output_name){
  df[output_name] <- df %>% select(all_of(electrode_list)) %>% rowMeans(na.rm = T)
  return(df)
}

mse_fix_outlier_avgTrials_wide <- combine_electrodes(df = mse_fix_outlier_avgTrials_wide, electrode_list = c("F7_Var1", "F8_Var1"), output_name = "vlpfc_Var1")
mse_fix_outlier_avgTrials_wide <- combine_electrodes(df = mse_fix_outlier_avgTrials_wide, electrode_list = c("AF5_Var1", "AF6_Var1", "F3_Var1", "F4_Var1", "F5_Var1", "F6_Var1"), output_name = "dlpfc_Var1")
mse_fix_outlier_avgTrials_wide <- combine_electrodes(df = mse_fix_outlier_avgTrials_wide, electrode_list = c("AF1_Var1", "AF2_Var1", "F1_Var1", "F2_Var1"), output_name = "spfc_Var1")
mse_fix_outlier_avgTrials_wide <- combine_electrodes(df = mse_fix_outlier_avgTrials_wide, electrode_list = c("FC1_Var1", "FC2_Var1", "FC3_Var1", "FC4_Var1"), output_name = "motor_Var1")

mse_fix_outlier_avgTrials_wide <- mse_fix_outlier_avgTrials_wide %>%
  mutate(across(c("vlpfc_Var1", "dlpfc_Var1", "spfc_Var1", "motor_Var1"), naoutlier)) %>% ungroup()


saveRDS(mse_fix_outlier_avgTrials_wide, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/mse_readyforMyelin_fix.RDS")

## Superficial and deep cortex R1 in areas of EEG electrodes ----
SGIGmyelin.EEGatlas.7T <- readRDS("/Volumes/Hera/Projects/corticalmyelin_development/BIDS/derivatives/surface_metrics/SGIGR1_electrodeatlas_finalsample.RDS") #generated with /surface_metrics/surface_measures/extract_depthdependent_R1.R

SGIGmyelin.EEGatlas.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- combine_electrodes(df = depth, electrode_list = c("F7", "F8"), output_name = "vlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF5", "AF6", "F3", "F4", "F5", "F6"), output_name = "dlpfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("AF1", "AF2", "F1", "F2"), output_name = "spfc_R1")
  depth <- combine_electrodes(df = depth, electrode_list = c("FC1", "FC2", "FC3", "FC4"), output_name = "motor_R1")
  depth <- depth %>% select(subject_id, session_id, subses, lunaid, visitno, eeg.date, age, sex, vlpfc_R1, dlpfc_R1, spfc_R1, motor_R1)
})

mse_fix_outlier_avgTrials_wide$lunaid <- as.integer(mse_fix_outlier_avgTrials_wide$lunaid)

myelin.mse.7T <- lapply(SGIGmyelin.EEGatlas.7T, function(depth){
  depth <- left_join(depth, mse_fix_outlier_avgTrials_wide, by = c("lunaid", "visitno"))
  return(depth)
})

saveRDS(myelin.mse.7T, "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/R1_mse_fix.RDS")





