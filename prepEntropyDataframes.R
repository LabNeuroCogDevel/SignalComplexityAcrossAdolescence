
# Load Libraries 

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)
library(eegUtils)
library(tvem)
library(interactions)
library(akima)


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) {ifelse(outliers(x), NA, x)}

# Load Dataframes ----

entropy6 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay6.csv') %>% select(c(-type)) %>% mutate(epoch = 'delayLength6')
entropy8 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay8.csv') %>% select(c(-type)) %>% mutate(epoch = 'delayLength8')
entropy10 <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay10.csv') %>% select(c(-type)) %>% mutate(epoch = 'delayLength10')
fix <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_fix.csv') %>% select(c(-type)) %>% mutate(epoch = 'Fix')

allDelayEntropy <- rbind(entropy6, entropy8) %>% rbind(., entropy10)%>% rbind(., fix)

# Entropy Outlier Detection ----
entropy_outlier <- allDelayEntropy %>% group_by(Subject, epoch) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% ungroup()

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

# Entropy averaged across all electrodes ----
entropyAgeAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, Var1, age) ~ Subject + visitno + epoch, data = entropyAge, FUN = mean)
write.csv(entropyAgeAvg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_MSE_delayLengths_allChansAvg.csv')

# Outlier detection on averaged electrodes ----
entropyAgeAvg_outlier <- entropyAgeAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,13,16,19,Inf), labels = c('10-12','13-15','16-18','Adults'))) %>% mutate(timeScale = as.numeric(timeScale))

write.csv(entropyAgeAvg_outlier_long, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_MSE_delayLengths_allChansAvg_outlierDect.csv')

## Find max entropy for every electrode ----
# Find each subjects max entropy value on the MSE curve 

entropyAge_long <- entropyAge %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))

subs = unique(entropyAge_long$Subject)
delayLengths = unique(entropyAge_long$epoch)
visits <- unique(na.omit(entropyAge_long$visitno))

maxValuesAllChans <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  
  subData <- entropyAge_long %>% filter(Subject == subject)
  chans = unique(subData$labels)
  
  for (delay in delayLengths) {
    
    subDelayData <- subData %>% filter(epoch == delay)
    
    for (chan in chans) {
      # Check if the subject has the visit number in the dataset
      
      subChanData <- subDelayData %>% filter(labels == chan)
      maxSubchanEntropy <- max(subChanData$MSx)
      
      subChanInfo <- data.frame(Subject = subject, labels = chan, epoch = delay, maxEntropy = maxSubchanEntropy)
      maxValuesAllChans <- rbind(maxValuesAllChans, subChanInfo)
    }
  }
}

maxValuesAllChansAge <- merge(maxValuesAllChans, ageValues, by = c('Subject')) %>% merge(., entropyAge %>% select("Subject", "labels", "Var1"), by = c("Subject", "labels"))

# Average Max Entropy ----

maxEntropyAvg <- aggregate(cbind(maxEntropy, Var1, age) ~ Subject + epoch, data = maxValuesAllChansAge, FUN = mean) 

# Outlier detection on averaged electrodes ----
maxEntropyAvg_outlier <- maxEntropyAvg %>% group_by(Subject) %>%
  mutate(across(c("maxEntropy", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))

write.csv(maxEntropyAvg_outlier, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_maxEntropy_delayLengths_allChansAvg.csv')
