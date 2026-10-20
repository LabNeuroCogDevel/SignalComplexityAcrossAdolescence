
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
library(extrafont)
extrafont::loadfonts()
library(conflicted)
library(eegkit)

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)



conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()


# PSD plot ----

psd_allElectrodes <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/attic/Results/PSDtable_allElectrodes.csv')

psd_allElectrodes <- psd_allElectrodes %>% 
  rename(urchan = Channel)

chanLocs <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Cortical_SNR_Development/resources/ChannelLocs.csv')%>% select(-type)

psd_allElectrodes <- merge(psd_allElectrodes, chanLocs %>% 
                             select('urchan', 'labels'), by = c('urchan'))

psd_allElectrodes <- psd_allElectrodes %>%
  mutate(region = case_when(
    startsWith(labels, "F") ~ "Frontal",
    startsWith(labels, "A") ~ "Frontal",
    startsWith(labels, "Cz") | startsWith(labels, "C1")  | startsWith(labels, "C2") ~ "Frontal",
    startsWith(labels, "C3") | startsWith(labels, "C4")  | startsWith(labels, "C5") | startsWith(labels, "C6") | startsWith(labels, "CPz") ~ "SensoryMotor",
    startsWith(labels, "CP") ~ "SensoryMotor",
    startsWith(labels, "P") ~ "Parietal",
    startsWith(labels, "O") ~ "Occipital",
    startsWith(labels, "I") ~ "Occipital",
    startsWith(labels, "T") ~ "Temporal"
  ))

psd_allElectrodes <- psd_allElectrodes %>% filter(region != "Temporal") %>% filter(region != "SensoryMotor")

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
subjectInfo <- merge7t[c("lunaid","eeg.date","visitno","eeg.age", "sex")]
subjectInfo$Subject <- paste(subjectInfo$lunaid, subjectInfo$eeg.date, sep = "_")
psd_allElectrodes <- merge(psd_allElectrodes, subjectInfo, 'Subject')


psd_allElectrodes_avgRegions <- psd_allElectrodes %>%
  mutate(ageGroup = cut(eeg.age, c(0,18,Inf), labels = c('Adol','Adults')))%>%
  group_by(freqs, region, ageGroup) %>%
  summarize(
    power = mean(power, na.rm = TRUE))


psd_long_filtered <- psd_allElectrodes %>% 
  mutate(ageGroup = cut(eeg.age, c(0,18,Inf), labels = c('Adol','Adults')))%>%
  filter(!is.infinite(power) & power > 0 & !is.na(power)) %>%
  group_by(ageGroup, region, freqs) %>% 
  mutate(z = scale(power)[,1],
         f = freqs,
         Band = cut(freqs, c(1,4,8,13,30,70), labels = c('Delta (1-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Gamma (30-48 Hz)'))
  ) %>%
  ungroup() %>%
  filter(abs(z) < 2)


psd_byGroup <- psd_long_filtered %>% 
  mutate(ageGroup = cut(eeg.age, c(0,18,Inf), labels = c('Adol','Adults')))%>%
  filter(!is.infinite(power) & power > 0 & !is.na(power)) %>%
  group_by(ageGroup, freqs, region) %>% 
  dplyr::summarize(meanPower = mean(power, na.rm=T), 
                   medPower = median(power, na.rm=T),
                   power_sd = sd(power, na.rm=T),
                   n = n()) %>%
  mutate(f = freqs,
         se = power_sd / sqrt(n))



lunaize(ggplot(data = psd_byGroup %>% filter(f>=4 & f<59) %>% arrange(f) %>% mutate(Age = ageGroup), 
                          aes(x = f, y = meanPower, color = region, linetype = Age)) +
                     geom_line() +
                     geom_ribbon(aes(ymin=meanPower-se, ymax=meanPower+se, fill=region), alpha=0.5, linewidth=0) +
                     scale_y_log10() +
                     labs(x = 'Frequency (Hz)', y = 'Power')) + scale_fill_manual(values=c("#012b0d", "#90b8a6", "#27763d"))+
  scale_color_manual(values=c("#012b0d", "#90b8a6", "#27763d"))+
  theme(text = element_text(size = 25),legend.position="right",
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = -20, l = 0)),
        axis.text.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = -25)), aspect.ratio = 1,axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -2),plot.margin = margin(0.5,0.5,0.5,0.5,"cm")) + facet_wrap(~region)



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/fooofFig/PSDs.pdf", device = "pdf", dpi = 500, width = 18, height = 6)



task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/taskRest.csv') %>%
  filter(epoch == 'restEyesOpen')

fooofAllChannels <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsAllChannelsFooofMeasures_20230911.csv')

fooofAllChannels <- fooofAllChannels %>%
  mutate(region = case_when(
    startsWith(Channel, "F") ~ "Frontal",
    startsWith(Channel, "A") ~ "Frontal",
    startsWith(Channel, "Cz") | startsWith(Channel, "C1")  | startsWith(Channel, "C2") ~ "Frontal",
    startsWith(Channel, "C3") | startsWith(Channel, "C4")  | startsWith(Channel, "C5") | startsWith(Channel, "C6") | startsWith(Channel, "CPz") ~ "SensoryMotor",
    startsWith(Channel, "CP") ~ "SensoryMotor",
    startsWith(Channel, "P") ~ "Parietal",
    startsWith(Channel, "O") ~ "Occipital",
    startsWith(Channel, "I") ~ "Occipital",
    startsWith(Channel, "T") ~ "Temporal"
  ))

fooofAllChannels <- fooofAllChannels %>% filter(region != "Temporal") %>% filter(region != "SensoryMotor")
                                                                                 
fooofAllChannels_avgRegions <- fooofAllChannels %>%
  group_by(Subject, age, sex, visitno, region, Condition) %>%
  summarize(
    Exponent = mean(Exponent, na.rm = TRUE),
    Offset = mean(Offset, na.rm = TRUE),
    age = mean(age, na.rm = TRUE)) %>% 
  mutate(ageGroup = cut(age, c(0,18,Inf), labels = c('Adol','Adults'))) %>% 
  separate(Subject, c('lunaid','vdate')) %>%
  filter(Condition == 'eyesOpen')


fooofAllChannels_avgRegions <- fooofAllChannels_avgRegions %>% group_by(lunaid) %>%
  mutate(across(c("Exponent","Offset"), naoutlier))

acw_fooof <- merge(task_rest, fooofAllChannels_avgRegions, by = c('lunaid', 'visitno', 'age','sex','region'))


lunaize(ggplot(data = acw_fooof, aes(x = ACW_50, y = Exponent, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  geom_smooth(aes(group=1), method='lm', alpha=0.5, linewidth = 1.5)+ facet_wrap(~region)+
  scale_color_manual(values = c("#012b0d", "#90b8a6", "#27763d")) + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/fooofFig/ACW_Exponent.pdf", device = "pdf", dpi = 500, width = 18, height = 6)


acw_fooof %>% 
  filter(region == 'Frontal') %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))


acw_fooof %>% 
  filter(region == 'Parietal') %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))

acw_fooof %>% 
  filter(region == 'Occipital') %>% 
  summarise(correlation = cor(ACW_50, Exponent, use = 'complete.obs'))


lm.model <- lmerTest::lmer(Exponent ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Frontal'))
summary(lm.model)

lm.model <- lmerTest::lmer(Exponent ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Parietal'))
summary(lm.model)

lm.model <- lmerTest::lmer(Exponent ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Occipital'))
summary(lm.model)



lunaize(ggplot(data = acw_fooof, aes(x = ACW_50, y = Offset, color=region))) + geom_point(alpha =0.3) +
  geom_line(aes(group= interaction(lunaid)), alpha = 0.2) +
  geom_smooth(aes(group=1), method='lm', alpha=0.5, linewidth = 1.5)+ facet_wrap(~region)+
  scale_color_manual(values = c("#012b0d", "#90b8a6", "#27763d")) + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = c("black")),
    axis.title.x = element_text(size = 20, color = c("black")),
    axis.title.y = element_text(size = 20, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white')) + facet_wrap(~region)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/fooofFig/ACW_Offset.pdf", device = "pdf", dpi = 500, width = 18, height = 6)


acw_fooof %>% 
  filter(region == 'Frontal') %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))

acw_fooof %>% 
  filter(region == 'Parietal') %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))

acw_fooof %>% 
  filter(region == 'Occipital') %>% 
  summarise(correlation = cor(ACW_50, Offset, use = 'complete.obs'))



lm.model <- lmerTest::lmer(Offset ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Frontal'))
summary(lm.model)

lm.model <- lmerTest::lmer(Offset ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Parietal'))
summary(lm.model)

lm.model <- lmerTest::lmer(Offset ~ ACW_50 + age + (1|lunaid), data = acw_fooof %>% filter(region == 'Occipital'))
summary(lm.model)
