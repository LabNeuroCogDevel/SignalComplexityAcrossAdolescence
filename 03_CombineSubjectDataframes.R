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

# Delay length 6 ----
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

# Delay length 8 ----

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

# Delay length 10 ----
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

# Fix ----
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


# Rest ----
## Two Second Epochs ----

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



