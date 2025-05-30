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

conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MonkeyData/samplingRate500/individual_files_delay/")

# List all CSV files in the directory
csv_files <- list.files()

# Initialize an empty data frame to store the combined data
combined_data_delay <- data.frame()

# Loop through each CSV file, read it, and combine it
# can also be done by readr::read_csv(Sys.glob("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MonkeyData/samplingRate500/individual_files_delay/*.csv"))

for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data_delay <- rbind(combined_data_delay, data)
}

combined_data_delay_outlier <- combined_data_delay %>% group_by(Monkey, BehavSess) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

combined_data_delay_outlier$epoch <- 'delay'



setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MonkeyData/samplingRate500/individual_files_fix/")

# List all CSV files in the directory
csv_files <- list.files()

# Initialize an empty data frame to store the combined data
combined_data_fix <- data.frame()

# Loop through each CSV file, read it, and combine it
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data_fix <- rbind(combined_data_fix, data)
}

combined_data_fix_outlier <- combined_data_fix %>% group_by(Monkey, BehavSess) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

combined_data_fix_outlier$epoch <- 'fix'

Delay_fix <- rbind(combined_data_fix_outlier, combined_data_delay_outlier)

Delay_fix$age_years <- Delay_fix$chronological_age / 365.25
Delay_fix$age_months <- Delay_fix$chronological_age / 30.44
Delay_fix$age_from_mature <- (Delay_fix$chronological_age) - (Delay_fix$MatureAge)
Delay_fix$months_from_mature <- Delay_fix$age_from_mature / 30.44


Delay_fix$Monkey <- as.factor(Delay_fix$Monkey)
Delay_fix$Channel <- as.factor(Delay_fix$Channel)
Delay_fix$BehavSess <- as.numeric(Delay_fix$BehavSess)
Delay_fix$Trial <- as.numeric(Delay_fix$Trial)
Delay_fix$NeuronArea <- as.factor(Delay_fix$NeuronArea)
Delay_fix$matureage_years <- Delay_fix$MatureAge / 365.25
Delay_fix$epoch <- factor(Delay_fix$epoch, levels = c('fix','delay'), ordered = TRUE)


my_palette <- c(
  "#012B0D",  # your very dark green
  "#0F4A23",  # added: deep green with blue undertone
  "#1C6530",  # added: forest green
  "#27763D",  # your mid green
  "#4D8D62",  # added: mossy green
  "#729D83",  # added: desaturated olive-green
  "#90B8A6"   # your muted light green
)

lunaize(ggplot(data = Delay_fix %>% filter(Monkey != 'QUA'), aes(x = age_years, y = ACW_50, color=Monkey, linetype=epoch))) +
  geom_smooth(aes(group=Monkey), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.1, linewidth = 1.5) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = F), alpha=0.1, linewidth = 2, linetype = 2, color='black') +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + 
  scale_color_manual(values = my_palette)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/monkeyFig/acwAcrossAge.pdf", device = "pdf", dpi = 500, width = 9, height = 6)


gam.model <-  mgcv::bam(ACW_50 ~ s(age_years, k = 3, by = Monkey) + s(BehavSess, k=3) + s(Trial, k =3) + s(Channel, by=Monkey, bs='re') + 
                          s(NeuronArea, bs='re') + s(Monkey, bs='re') + s(age_years, k = 3, by = epoch) + s(epoch, bs='re'), 
                        data = Delay_fix %>% filter(Monkey != 'QUA'), discrete = T)
summary(gam.model)


# bam.model <- mgcv::bam(ACW_50 ~ 
#   s(age_years, k = 3) +
#   s(age_years, k=3, by=epoch) +
#   BehavSess + #controls for differences in ACW by session.
#   Channel + #controls for differences in ACW by channel.
#   NeuronArea + #controls for differences in ACW by brain area.
#   epoch + #controls for difference in ACW by epoch.
#   s(age_years, k = 3, by = Monkey) + #Allows each monkey to have own age smooth
#   s(Trial, k = 3) #Controls for difference in ACW across trials.  
# + s(Monkey, bs = "re"), #Random part. Each monkey has own "starting point" for ACW
# data = Delay_fix %>% filter(Monkey != 'QUA'), discrete = T)
# 
# summary(bam.model)



predictDF <- ggeffects::ggpredict(gam.model, terms=c('epoch'))
ggeffects::test_predictions(predictDF, terms='epoch')


plot(ggeffects::ggpredict(gam.model, terms=c('age_years', 'Monkey')))
plot(ggeffects::ggpredict(gam.model, terms=c('Channel', 'Monkey')))
plot(ggeffects::ggpredict(gam.model, terms=c('NeuronArea')))
plot(ggeffects::ggpredict(gam.model, terms=c('BehavSess')))
plot(ggeffects::ggpredict(gam.model, terms=c('age_years', 'matureage_years_diff')))
plot(ggeffects::ggpredict(gam.model, terms=c('matureage_years_diff', 'Monkey')))
plot(ggeffects::ggpredict(gam.model, terms=c('epoch', 'age_years')))



lunaize(ggplot(data = Delay_fix %>% filter(Monkey != 'QUA'), aes(x = epoch, y = ACW_50, fill = interaction(factor(epoch, levels = c("delay", "fix"))))) +
          geom_boxplot(alpha = 0.8, color = "black", width = 0.5, outlier.shape = NA) +
          labs(x = "Epoch", y = "ACW", fill = "Task") +
          scale_fill_manual(values = c("#012b0d", "#27763d", "#012b0d", "#27763d", "#012b0d", "#27763d"))) + coord_cartesian(ylim = c(0, 0.025)) 

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/monkeyFig/delayVsFix.pdf", device = "pdf", dpi = 500, width = 9, height = 6)


Delay_fix_wide <- Delay_fix %>%
  select(Monkey, BehavSess, Channel, Trial, epoch, ACW_50) %>%
  pivot_wider(
    names_from = epoch,
    values_from = ACW_50,
    names_prefix = "ACW_50_"
  )

t.test(Delay_fix_wide$ACW_50_delay, Delay_fix_wide$ACW_50_fix, paired = T)

Delay_fix_wide %>%
  group_by(Monkey, BehavSess) %>%
  summarise(t_test = list(t.test(ACW_50_delay, ACW_50_fix, paired = TRUE))) %>%
  mutate(
    p_value = purrr::map_dbl(t_test, ~ .x$p.value),
    t_stat  = purrr::map_dbl(t_test, ~ .x$statistic),
    df      = purrr::map_dbl(t_test, ~ .x$parameter)
  )

lunaize(ggplot(data = Delay_fix %>% filter(Monkey != 'QUA'), aes(x = Monkey, y = ACW_50, fill = interaction(factor(epoch, levels = c("fix", "delay"))))) +
          geom_boxplot(alpha = 0.8, color = "black", width = 0.5, outlier.shape = NA) +
          labs(x = "Epoch", y = "ACW", fill = "Task") +
          scale_fill_manual(values = c("#012b0d", "#27763d", "#012b0d", "#27763d", "#012b0d", "#27763d"))) + coord_cartesian(ylim = c(0, 0.025)) 

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/monkeyFig/delayVsFix_eachMonkey.pdf", device = "pdf", dpi = 500, width = 9, height = 6)








# Average across trials ----
combined_data_outlier_avgTrials <- combined_data_outlier %>%
  group_by(Monkey, BehavSess, Channel, MatureAge, chronological_age, group, mature_group, age_group, age_years, months_from_mature, NeuronArea) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE))


combined_data_outlier_avgTrials <- combined_data_outlier_avgTrials %>% group_by(Monkey) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

lunaize(ggplot(data = combined_data_outlier_avgTrials%>% filter(Monkey != 'QUA'), aes(x = age_years, y = ACW_50, color=Monkey)))  +
  geom_smooth(aes(group=Monkey), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.1, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))


gam.model <-  mgcv::bam(ACW_50 ~ s(age_years, k = 3, by = Monkey) + BehavSess + Channel + Monkey, data = combined_data_outlier_avgTrials %>% filter(Monkey != 'QUA'), discrete = T)
summary(gam.model)

plot(ggeffects::ggpredict(gam.model, terms=c('age_years', 'Monkey')))


# Average across behavioral sessions ----

combined_data_outlier_avgTrials_avgChans <- combined_data_outlier_avgTrials %>%
  group_by(Monkey, BehavSess, MatureAge, chronological_age, group, mature_group, age_group, age_years, months_from_mature, NeuronArea) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE)) %>% group_by() %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

lunaize(ggplot(data = combined_data_outlier_avgTrials_avgChans, aes(x = age_years, y = ACW_50, color=Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) 


gam.model <-  mgcv::gamm(ACW_50 ~ s(age_years, k = 3) , data = combined_data_outlier_avgTrials_avgChans, random=list(Monkey=~1))
summary(gam.model$gam)


lunaize(ggplot(data = combined_data_outlier_avgTrials_avgBehav, aes(x = age_years, y = ACW_0, color=Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))


gam.model <-  mgcv::gamm(ACW_50 ~ s(age_years, k = 3) , data = combined_data_outlier_avgTrials_avgBehav, random=list(Monkey=~1))
summary(gam.model$gam)

# Average across channels ----

combined_data_outlier_avgTrials_avgChan <- combined_data_outlier_avgTrials %>%
  group_by(Monkey, BehavSess, MatureAge, chronological_age, group, mature_group, age_group, age_years, months_from_mature) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE)) %>% group_by(Monkey) %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()


lunaize(ggplot(data = combined_data_outlier_avgTrials_avgChan, aes(x = age_years, y = ACW_0, color = Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) 

gam.model <-  mgcv::gamm(ACW_50 ~ s(months_from_mature, k = 3) , data = combined_data_outlier_avgTrials_avgBehav_avgChan, random=list(Monkey=~1))
summary(gam.model$gam)


lunaize(ggplot(data = combined_data_outlier_avgTrials_avgBehav_avgChan, aes(x = age_years, y = ACW_50, color = Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))


gam.model <-  mgcv::gamm(ACW_50 ~ s(age_years, k = 3) , data = combined_data_outlier_avgTrials_avgBehav_avgChan, random=list(Monkey=~1))
summary(gam.model$gam)


# Avg channels within a behavioral ----
combined_data_outlier_avgTrials_avgchans <- combined_data_outlier_avgTrials %>%
  group_by(Monkey, BehavSess, MatureAge, chronological_age, group, mature_group, age_group, age_years, months_from_mature, NeuronArea) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE)) %>% group_by() %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

lunaize(ggplot(data = combined_data_outlier_avgTrials_avgchans, aes(x = age_years, y = ACW_50, color = Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))


# avg behav after channels ----

combined_data_outlier_avgTrials_avgchans_avgBehav <- combined_data_outlier_avgTrials_avgchans %>%
  group_by(Monkey, MatureAge, chronological_age, group, mature_group, age_group, age_years, months_from_mature, NeuronArea) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE)) %>% group_by() %>%
  mutate(across(c("ACW_50"), naoutlier)) %>% ungroup()

lunaize(ggplot(data = combined_data_outlier_avgTrials_avgchans_avgBehav, aes(x = age_years, y = ACW_50, color = Monkey))) + geom_point(alpha =0.3) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))


# Define half-year bins
breaks <- seq(-31, 20, by = 5)
labels <- paste0(head(breaks, -1), "–", tail(breaks, -1))

# Add a new column with the grouped ages
combined_data_outlier_avgTrials_avgBehav_avgChan$age_month_group <- cut(
  combined_data_outlier_avgTrials_avgBehav_avgChan$months_from_mature,
  breaks = breaks,
  labels = labels,
  right = FALSE
)

acw_groupedbyage <- combined_data_outlier_avgTrials_avgBehav_avgChan %>%
  group_by(Monkey, age_month_group) %>%
  summarize(
    ACW_0 = mean(ACW_0, na.rm = TRUE),
    ACW_50 = mean(ACW_50, na.rm = TRUE))

lunaize(ggplot(data = acw_groupedbyage, aes(x = age_month_group, y = ACW_50))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(Monkey)), alpha = 0.2) +
  geom_smooth(aes(group=1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.3, linewidth = 1.5)  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'))
