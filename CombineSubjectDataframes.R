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

# Initialize directory and empty dataframes ----
chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/ChannelLocs.csv')

# Set your working directory to the folder containing your CSV files
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/individual_subject_files/")


# delay length 6 ----

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay6.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, chanLocs)
  
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay6.csv', row.names = F)


# delay length 8 ----

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay8.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, chanLocs)
  
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay8.csv', row.names = F)


# delay length 10 ----

# List all CSV files in the directory
csv_files <- list.files(pattern = "delay10.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  subdata <- data[1:64, ] 
  subdata <- cbind(subdata, chanLocs)
  
  combined_data <- rbind(combined_data, subdata)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_delay10.csv', row.names = F)


# fix ----

# List all CSV files in the directory
csv_files <- list.files(pattern = "fix.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  subdata <- data[1:64, ] 
  subdata <- cbind(subdata, chanLocs)
  
  combined_data <- rbind(combined_data, subdata)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_allChans_MSE20_fix.csv', row.names = F)

