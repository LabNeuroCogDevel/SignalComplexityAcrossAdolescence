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

## ACW ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/individual_subject_files/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "ACW.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/allSubjects_ACW.csv', row.names = F)



## ROIs ----

datapath <- "/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/"
outpath <- "/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/"

all_files <- list.files(path = datapath, full.names = TRUE, recursive = TRUE)
setfiles0 <- grep("P.*/[Rr]est/channelROIs.csv", all_files, value = TRUE)

# Store file paths in a list (similar to cell array in MATLAB)
setfiles <- as.list(setfiles0)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in setfiles) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, file = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/allSubjectsROIs.csv', row.names = F)
