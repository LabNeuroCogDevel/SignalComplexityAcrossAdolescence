# Libraries ----

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

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# Initialize directory and empty dataframes ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$Subject <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")

# delay length 6 ----
## MSE ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay6.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}
combined_data <- combined_data %>% merge(., ageValues, by = "Subject")

combined_data_outlier <- combined_data %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10","MSx11", "MSx12", "MSx13", "MSx14", "MSx15", "MSx16", 
                  "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% ungroup()

combined_data_outlier$secondsEpoch <- as.factor(combined_data_outlier$secondsEpoch)
combined_data_outlier$Subject <- as.factor(combined_data_outlier$Subject)


combined_data_outlier_avg <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(combined_data_outlier_avg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgChans_MSE20_delay6.csv', row.names = F)

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')

combined_data_outlier_avgTrials <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

combined_data_outlier_avgTrials_chans <- merge(chanLocs, combined_data_outlier_avgTrials, by ="channel")

write.csv(combined_data_outlier_avgTrials_chans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay6.csv', row.names = F)



## Power Spectrum ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/MGS/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay6.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/allSubjects_allChans_powerSpec_delay6.csv', row.names = F)


## ACW ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*delay6.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay6.csv', row.names = F)

# delay length 8 ----
## MSE ----
# List all CSV files in the directory
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/individual_subject_files/")

csv_files <- list.files(pattern = "delay8.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

combined_data <- combined_data %>% merge(., ageValues, by = "Subject")

combined_data_outlier <- combined_data %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10","MSx11", "MSx12", "MSx13", "MSx14", "MSx15", "MSx16", 
                  "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% ungroup()

combined_data_outlier$secondsEpoch <- as.factor(combined_data_outlier$secondsEpoch)
combined_data_outlier$Subject <- as.factor(combined_data_outlier$Subject)


combined_data_outlier_avg <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(combined_data_outlier_avg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/allSubjects_avgChans_MSE20_delay8.csv', row.names = F)


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')

combined_data_outlier_avgTrials <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

combined_data_outlier_avgTrials_chans <- merge(chanLocs, combined_data_outlier_avgTrials, by ="channel")
write.csv(combined_data_outlier_avgTrials_chans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay8.csv', row.names = F)


## Power Spectrum ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/MGS/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay8.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed

  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/allSubjects_allChans_powerSpec_delay8.csv', row.names = F)

## ACW ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*delay8.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay8.csv', row.names = F)

# delay length 10 ----
## MSE ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay10.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

combined_data <- combined_data %>% merge(., ageValues, by = "Subject")

combined_data_outlier <- combined_data %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10","MSx11", "MSx12", "MSx13", "MSx14", "MSx15", "MSx16", 
                  "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% ungroup()

combined_data_outlier$secondsEpoch <- as.factor(combined_data_outlier$secondsEpoch)
combined_data_outlier$Subject <- as.factor(combined_data_outlier$Subject)


combined_data_outlier_avg <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

write.csv(combined_data_outlier_avg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/allSubjects_avgChans_MSE20_delay10.csv', row.names = F)


chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')

combined_data_outlier_avgTrials <- combined_data_outlier %>%
  group_by(Subject, secondsEpoch, channel, age, sex, visitno) %>%
  summarize(
    MSx1 = mean(MSx1, na.rm = TRUE),
    MSx2 = mean(MSx2, na.rm = TRUE),
    MSx3 = mean(MSx3, na.rm = TRUE),
    MSx4 = mean(MSx4, na.rm = TRUE),
    MSx5 = mean(MSx5, na.rm = TRUE),
    MSx6 = mean(MSx6, na.rm = TRUE),
    MSx7 = mean(MSx7, na.rm = TRUE),
    MSx8 = mean(MSx8, na.rm = TRUE),
    MSx9 = mean(MSx9, na.rm = TRUE),
    MSx10 = mean(MSx10, na.rm = TRUE),
    MSx11 = mean(MSx11, na.rm = TRUE),
    MSx12 = mean(MSx12, na.rm = TRUE),
    MSx13 = mean(MSx13, na.rm = TRUE),
    MSx14 = mean(MSx14, na.rm = TRUE),
    MSx15 = mean(MSx15, na.rm = TRUE),
    MSx16 = mean(MSx16, na.rm = TRUE),
    MSx17 = mean(MSx17, na.rm = TRUE),
    MSx18 = mean(MSx18, na.rm = TRUE),
    MSx19 = mean(MSx19, na.rm = TRUE),
    MSx20 = mean(MSx20, na.rm = TRUE),
    Var1 = mean(Var1, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% separate(Subject, c('lunaid','vdate')) %>% mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults')))

names(chanLocs)[names(chanLocs) == "urchan"] <- "channel"

combined_data_outlier_avgTrials_chans <- merge(chanLocs, combined_data_outlier_avgTrials, by ="channel")
write.csv(combined_data_outlier_avgTrials_chans, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/twoSecondEpochs/allSubjects_avgTrials_MSE20_delay10.csv', row.names = F)


## Power Spectrum ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/MGS/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay10.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed

  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/allSubjects_allChans_powerSpec_delay10.csv', row.names = F)


## ACW ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*delay10.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/twoSecondSegments/allSubjects_allChans_acw_delay10.csv', row.names = F)

# fix ----
## MSE ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "fix0.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/MGS/allSubjects_allChans_mse_fix.csv', row.names = F)

## Power Spectrum ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/MGS/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "fix.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/allSubjects_allChans_powerSpec_fix.csv', row.names = F)


## ACW ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*fix0.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MGS/allSubjects_allChans_acw_fix.csv', row.names = F)


# rest ----
## MSE ----
### eyes closed ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/individual_subject_files/timeScale20/eyesClosed/")


# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesClosed.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  subdata <- data[1:64, ] 
  subdata <- cbind(subdata, chanLocs)
  
  combined_data <- rbind(combined_data, subdata)
}
combined_data$epoch <- "restEyesClosed"

### eyes open ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/entropy/individual_subject_files/timeScale20/eyesOpen/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "Entropy.csv")

# Initialize an empty data frame to store the combined data
combined_data_open <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  subdata <- data[1:64, ] 
  subdata <- cbind(subdata, chanLocs)
  
  combined_data_open <- rbind(combined_data_open, subdata)
}
combined_data_open$epoch <- "restEyesOpen"

open_closed <- rbind(combined_data, combined_data_open)
  
write.csv(open_closed, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/entropy/allSubjects_allChans_MSE20_rest.csv', row.names = F)

## ACW ----
### eyes closed ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*eyesClosed.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesClosed.csv', row.names = F)


### eyes open ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "acw.*eyesOpen.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/allSubjects_allChans_acw_eyesOpen.csv', row.names = F)


## Power Spectrum ----
### eyes closed ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/Rest/")


# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesClosed.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}
combined_data$epoch <- "restEyesClosed"

### eyes open ----

# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesOpen.csv")

# Initialize an empty data frame to store the combined data
combined_data_open <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data_open <- rbind(combined_data_open, data)
}

combined_data_open$epoch <- "restEyesOpen"

open_closed <- rbind(combined_data, combined_data_open)

write.csv(open_closed, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/powerSpectrum/allSubjects_allChans_powerSpec_rest.csv', row.names = F)


## Two Second Epochs ----
### MSE ----
#### eyes closed ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/individual_subject_files/")


# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesClosed0.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}
combined_data$epoch <- "restEyesClosed"

#### eyes open ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesOpen0.csv")

# Initialize an empty data frame to store the combined data
combined_data_open <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data_open <- rbind(combined_data_open, data)
}
combined_data_open$epoch <- "restEyesOpen"

open_closed <- rbind(combined_data, combined_data_open)

write.csv(open_closed, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/mse/Rest/twoSecondSegments/allSubjects_allChans_MSE20_2seconds_rest.csv', row.names = F)



### ACW ----
#### eyes closed ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/individual_subject_files/")


# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesClosed0.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}
combined_data$epoch <- "restEyesClosed"

#### eyes open ----

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "eyesOpen0.csv")

# Initialize an empty data frame to store the combined data
combined_data_open <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data_open <- rbind(combined_data_open, data)
}
combined_data_open$epoch <- "restEyesOpen"

open_closed <- rbind(combined_data, combined_data_open)

write.csv(open_closed, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/Rest/twoSecondSegments/allSubjects_allChans_acw_2seconds_rest.csv', row.names = F)


