

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

# Load dataframes ----

entropyAgeAvg_outlier_long <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_MSE_delayLengths_allChansAvg_outlierDect.csv')
maxEntropyAvg_outlier <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/MGS_Entropy/allSubjects_maxEntropy_delayLengths_allChansAvg.csv')

# Plot entropy across timescales for each delay length ----

lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = timeScale, y = MSx, color = as.factor(epoch))) + geom_point() + 
          geom_smooth(aes(group = as.factor(epoch)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

# Plot AUC across age for each delay length ----

lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = age, y = Var1, color = as.factor(epoch))) + geom_point() + 
          geom_line(aes(group= interaction(lunaID, epoch)), alpha = 0.2) +
          geom_smooth(aes(group = as.factor(epoch)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


# Plot max entropy across age for each delay length ----

lunaize(ggplot(data = maxEntropyAvg_outlier, aes(x = age, y = maxEntropy, color = as.factor(epoch))) + geom_point() + 
          geom_line(aes(group= interaction(lunaID, epoch)), alpha = 0.2) +
          geom_smooth(aes(group = as.factor(epoch)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))



