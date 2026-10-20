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
ageValues <- merge7t[c("lunaid","rest.date","visitno","eeg.age", "sex")]
colnames(ageValues) <- c("lunaid", "visitDate","visitno", "age", "sex")
ageValues$SubjectID <- paste(ageValues$lunaid, ageValues$visitDate, sep = "_")

# ----
## ACW ----
### rest ----
setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/individual_subject_files_beamforming//")
regionOrder <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasserRegionOrder.csv')

# List all CSV files in the directory
csv_files <- list.files(pattern = "*rest.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}

colnames(combined_data)[colnames(combined_data) == "SubjectIDTable"] <- "SubjectID"

age_combined_data <- merge(combined_data, ageValues, by = 'SubjectID')

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

age_combined_data <- age_combined_data %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

age_combined_data <- merge(age_combined_data, glasser_lobes, by = 'orig_parcelname')


  
write.csv(age_combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beam_rest_epochs.csv', row.names = F)

age_combined_data <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beam_rest_epochs.csv')

glasserINTsRest <- age_combined_data %>% group_by(SubjectID, lobe) %>%
  mutate(across(c("ACW_50"), naoutlier))

glasserINTs_avgLobe <-  age_combined_data %>%  group_by(SubjectID, lobe, visitno, age) %>%
  summarise(INT_glasser = mean(ACW_50, na.rm = TRUE)) 

glasserINTs_avgLobe <- glasserINTs_avgLobe %>%
  mutate(region = recode(lobe,
                         "F" = "Frontal",
                         "O" = "Occipital",
                         "P" = "Parietal")) %>% separate(SubjectID, c('lunaid','vdate'))

lunaize(ggplot(data = glasserINTsRest, aes(x = age, y = ACW_50, color=lobe))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5)+ facet_wrap(~lobe)+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))



### Task activity ----
regionOrder <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasserRegionOrder.csv')

setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/individual_subject_files_beamforming/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "*delay6..csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '6';
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}



csv_files <- list.files(pattern = "*delay8..csv")

# Initialize an empty data frame to store the combined data
combined_data_8 <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '8';
  data <- cbind(data, regionOrder)
  combined_data_8 <- rbind(combined_data_8, data)
}





csv_files <- list.files(pattern = "*delay10.csv")

# Initialize an empty data frame to store the combined data
combined_data_10 <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '10';
  data <- cbind(data, regionOrder)
  combined_data_10 <- rbind(combined_data_10, data)
}


combine_all <- rbind(combined_data, combined_data_8) %>% rbind(., combined_data_10) %>% merge(., ageValues, by='SubjectID')

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

combine_all <- combine_all %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

combine_all <- merge(combine_all, glasser_lobes, by = 'orig_parcelname')


write.csv(combine_all, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_delayPeriods.csv', row.names = F)
 

### fix ----
# List all CSV files in the directory
csv_files <- list.files(pattern = "*fix.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

combine_all <- combined_data %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

combine_all <- merge(combine_all, glasser_lobes, by = 'orig_parcelname')

names(combine_all)[names(combine_all) == "SubjectIDTable"] <- "SubjectID"

age_combined_data <- merge(combine_all, ageValues, by = 'SubjectID')


write.csv(age_combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_fix.csv', row.names = F)





# ACW Varibaility ----
# List all CSV files in the directory
## rest ----
csv_files <- list.files(pattern = "*rest_variability.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}

colnames(combined_data)[colnames(combined_data) == "SubjectIDTable"] <- "SubjectID"

age_combined_data <- merge(combined_data, ageValues, by = 'SubjectID')

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

age_combined_data <- age_combined_data %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

age_combined_data <- merge(age_combined_data, glasser_lobes, by = 'orig_parcelname')



write.csv(age_combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beam_rest_epochs_variability.csv', row.names = F)


## delay ----
regionOrder <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasserRegionOrder.csv')

setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/individual_subject_files_beamforming/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "*delay6_variability.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '6';
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}



csv_files <- list.files(pattern = "*delay8_variability.csv")

# Initialize an empty data frame to store the combined data
combined_data_8 <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '8';
  data <- cbind(data, regionOrder)
  combined_data_8 <- rbind(combined_data_8, data)
}





csv_files <- list.files(pattern = "*delay10_variability.csv")

# Initialize an empty data frame to store the combined data
combined_data_10 <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data$delayLength = '10';
  data <- cbind(data, regionOrder)
  combined_data_10 <- rbind(combined_data_10, data)
}


combine_all <- rbind(combined_data, combined_data_8) %>% rbind(., combined_data_10) %>% mutate(SubjectID = SubjectIDTable) %>%
  merge(., ageValues, by='SubjectID') %>% select(-SubjectIDTable)

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

combine_all <- combine_all %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

combine_all <- merge(combine_all, glasser_lobes, by = 'orig_parcelname')


write.csv(combine_all, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_delayPeriods_variability.csv', row.names = F)



## fix ----
# List all CSV files in the directory
csv_files <- list.files(pattern = "*fix_variability.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  data <- cbind(data, regionOrder)
  combined_data <- rbind(combined_data, data)
}

glasser_lobes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser_lobes.csv')

combine_all <- combined_data %>%
  mutate(
    orig_parcelname = sub("\\s+[LR]$", "", OrigParcelName)
  ) 

combine_all <- merge(combine_all, glasser_lobes, by = 'orig_parcelname')

names(combine_all)[names(combine_all) == "SubjectIDTable"] <- "SubjectID"

age_combined_data <- merge(combine_all, ageValues, by = 'SubjectID')


write.csv(age_combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_fix_variability.csv', row.names = F)
