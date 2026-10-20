
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
library(kableExtra)
library(ggnewscale)
library(ggseg)
library(ggsegGlasser)


conflict_prefer("filter", "dplyr")  # Always use dplyr::filter()
conflicted::conflicts_prefer(lme4::lmer)

source("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code/behavioralStats.R")

lm.model.function <- function(input.df, region, dep_var, id_var){
  
  form <- as.formula(
    paste(dep_var, "~", id_var, "+ age + sex + epoch + (1|lunaid)")
  )
  
  lm.model <- lmerTest::lmer(
    form,
    data = input.df %>%
      filter(orig_parcelname == region)
  )
  
  lm.stats <- as.data.frame(summary(lm.model)$coefficients)
  lm.stats$term <- rownames(lm.stats)
  rownames(lm.stats) <- NULL
  
  return(lm.stats)
}


lm.interaction.model.function <- function(input.df, region, dep_var, id_var){
  
  form <- as.formula(
    paste(dep_var, "~", id_var, "*age + age + sex + epoch + (1|lunaid)")
  )
  
  lm.model <- lmerTest::lmer(
    form,
    data = input.df %>%
      filter(orig_parcelname == region)
  )
  
  lm.stats <- as.data.frame(summary(lm.model)$coefficients)
  lm.stats$term <- rownames(lm.stats)
  rownames(lm.stats) <- NULL
  return(lm.stats)
  
}


# Load in data ----

merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
behav <- merge7t[c("lunaid","eeg.date", "visitno","eeg.age", "eeg.vgsLatency_DelayAll","eeg.BestError_DelayAll", "eeg.BestError_sd_DelayAll","eeg.mgsLatency_DelayAll", "eeg.mgsLatency_sd_DelayAll","cantab.ssp_max.span","cantab.ssp_nerrors", "cantab.ssp_ntrials")]
colnames(behav) <- c("lunaid","visitdate", "visitno","age","vgsLatency","absBestError","absBestError_sd","mgsLatency","mgsLatency_sd", "SSP_maxSpan", "SSP_nErrors", "SSP_nTrials")

task_rest <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/allSubjects_glasserRegions_beamforming_taskRest_allregions_2.csv')

ACW_behav <- merge(task_rest, behav, by = c("lunaid", "visitno", "age")) 


# Absolute Best Error ----
regions <- unique(ACW_behav$orig_parcelname)
ACW_behav$lunaid <- as.factor(ACW_behav$lunaid)
ACW_behav$epoch <- as.factor(ACW_behav$epoch)

# Initialize lists to store results
maineffect_list <- list()
interactions_list <- list()

# Loop through each combination of epoch and region
  for (loc in regions) {
    
    lm.stats <- lm.model.function(input.df = ACW_behav, region = loc, dep_var = 'absBestError', id_var = 'ACW_50')
    
    # Create a data frame of stats and add epoch  info
   
    lm.stats.df <- data.frame(ACWestimate = lm.stats[2,1], ACWpValue = lm.stats[2,5], sexEstimate = lm.stats[4,1], sexpValue= lm.stats[4,5], epochEstimate = lm.stats[5,1], epochpValue = lm.stats[5,5]) %>% 
      mutate(orig_parcelname = loc)
      
    # Store results
    maineffect_list[[loc]] <- lm.stats.df
    
    
    lmInteraction.stats <- lm.interaction.model.function(input.df = ACW_behav, region = loc, dep_var = 'absBestError', id_var = 'ACW_50')
    
    lmInteraction.stats.df <- data.frame(estimate = lmInteraction.stats[6,1], pValue = lmInteraction.stats[6,5])%>% 
      mutate(orig_parcelname = loc)
  
    
    interactions_list[[loc]] <- lmInteraction.stats.df
    
  }


saveRDS(maineffect_list, file = "acw_behav_maineffects_linear_absBestError.rds")
saveRDS(interactions_list, file = "acw_behav_interactioneffects_linear_absBestError.rds")


maineffect_list <- readRDS("acw_behav_maineffects_linear_absBestError.rds")
interactions_list <- readRDS("acw_behav_interactioneffects_linear_absBestError.rds")


lmstats_combined <- bind_rows(maineffect_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

interaction_combined <- bind_rows(interactions_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )



glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

## run when remote ------------------------------------------------------------
# glasser <- glasser()
# glasser$data$hemi <- factor(glasser$data$hemi, levels = c('right','left'))
# glasser$data$side <- factor(glasser$data$side, levels = c('medial','lateral'))


lmstats_combined <- lmstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 

interaction_combined <- interaction_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

## run when remote ------------------------------------------------------------
# lmstats_combined <- lmstats_combined %>%
#   mutate(label = paste0(label, "_ROI"))
# 
# interaction_combined <- interaction_combined %>%
#   mutate(label = paste0(label, "_ROI"))


dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#f0faf3"   # very pale green
)


lmstats_combined <- lmstats_combined %>%
  mutate(
    correctedP = p.adjust(ACWpValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate)
  ) +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




interaction_combined <- interaction_combined %>%
  mutate(correctedP = p.adjust(pValue, method = "fdr")) %>%
  mutate(significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate)
  ) +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  (dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/accuracyInteraction.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 3.5, height = 3.5)

sigRegions <- interaction_combined %>% filter(significant == TRUE) 

sigRegionsNames <- sigRegions$orig_parcelname


sigRegionsAvg <- ACW_behav %>%
  filter(orig_parcelname %in% sigRegionsNames) %>%
  group_by(lunaid, visitno, epoch, age, sex, absBestError) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    .groups = "drop"
  )

lm.model <- lmerTest::lmer(
  "absBestError ~ ACW_50*age + sex + epoch + (1|lunaid)",
  data = sigRegionsAvg
)

plot(ggpredict(lm.model, terms=c('ACW_50', 'age')))
library(interactions)

jn <- interactions::johnson_neyman(
  lm.model,
  pred = ACW_50,
  modx = age
)

jn$plot +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  labs(
    title = NULL)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Supp/JohnsonNeymanAllsigregions.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 3.5, height = 3.5)


summary(lm.model)

plot(ggpredict(lm.model, terms=c('ACW_50', 'age')))

plotdata <- ggpredict(lm.model, terms = c("ACW_50", "epoch"), interval = "confidence")

plot(ggpredict(lm.model, terms = c("age" ,"ACW_50"), interval = "confidence"), show_data = T)



#### Neurosynth on all sig regions avg ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) 

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#001a08") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))

## Neurosyth correlation on estimate ----

interaction_combined_neurosynth <- merge(neurosynth.terms, interaction_combined, by='label')  %>% filter(significant == TRUE) 

term_cols <- setdiff(names(neurosynth.terms), "label")

# Correlate each term with estimate
cor_results <- lapply(term_cols, function(term) {
  
  x <- interaction_combined_neurosynth[[term]]
  y <- interaction_combined_neurosynth$estimate
  
  valid <- complete.cases(x, y)
  
  test <- cor.test(
    x[valid],
    y[valid],
    method = "pearson"
  )
  
  data.frame(
    term = term,
    r = unname(test$estimate),
    p_value = test$p.value,
    n = sum(valid)
  )
}) %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adjusted)

cor_results




###  p47r ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_p47r_ROI', 'L_p47r_ROI')) %>%
  group_by(lunaid, visitno, epoch) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  )

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#001a08"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#001a08")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 8, color = "#001a08"),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_47.pdf", device = "pdf", dpi = 500, width = 2, height = 2.5)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "R_p47r_ROI" = "#001a08",
    "L_p47r_ROI" = "#001a08"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_47.pdf", device = "pdf", dpi = 500, width = 18, height = 5)


###  dlpfc ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_a9-46v_ROI', 'L_a9-46v_ROI')) %>%
  group_by(lunaid, visitno, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  )%>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#001a08"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#001a08")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_dlpfc.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2)

lm.model <- lmerTest::lmer("absBestError ~ ACW_50*age + epoch + (1|lunaid)", data = ACW_behav_avg)

summary(lm.model)

pred <- (ggpredict(lm.model, terms=c('ACW_50[all]', 'age[14, 20, 26]')))



lunaize(
  ggplot() +
    # Raw data points
    geom_point(
      data = ACW_behav_avg,
      aes(x = ACW_50, y = absBestError, shape=epoch),
      color = "#001a08",
      alpha = 0.25,
      size = 0.6
    ) + 
    geom_line(data = ACW_behav_avg,
              aes(x = ACW_50, y = absBestError, shape=epoch, group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
    # Model-predicted confidence intervals
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = conf.low, ymax = conf.high,
          fill = group),
      alpha = 0.15,
      color = NA
    ) +
    # Interaction lines
    geom_line(
      data = pred,
      aes(x = x, y = predicted, color = group),
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      name = "Age (years)"
    ) +
    scale_fill_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      guide = "none"
    )) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      axis.line = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2),
      panel.grid = element_blank()
    ) +
    labs(
      x = "ACW (s)",
      y = "Accuracy (degrees)",
      color = "Age (years)"
    )



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_dlpfc_interaction.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)





library(interactions)

jn <- interactions::johnson_neyman(
  lm.model,
  pred = ACW_50,
  modx = age
)

jn$plot +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  labs(
    title = NULL)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_dlpfc_johnson.pdf", device = "pdf", dpi = 500, width = 1.7, height = 1.9)



ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "R_a9-46v_ROI" = "#001a08",
    "L_a9-46v_ROI" = "#001a08"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_47.pdf", device = "pdf", dpi = 500, width = 18, height = 5)



#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('R_a9-46v_ROI', 'L_a9-46v_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#001a08") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_dlpfc.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 2.1, height = 1.4)




###  motor ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI')) %>%
  group_by(lunaid, visitno, epoch) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#014421"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#014421")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 8, color = "#014421"),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_motor.pdf", device = "pdf", dpi = 500, width = 2, height = 2.5)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "R_4_ROI" = "#014421",
    "L_4_ROI" = "#014421"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_motor.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('R_4_ROI', 'L_4_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#001a08") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_motor.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 2.1, height = 1.4)




###  PF  ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_PF_ROI', 'L_PF_ROI')) %>%
  group_by(lunaid, visitno, epoch) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#014421"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#014421")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 8, color = "#014421"),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_PF_ROI.pdf", device = "pdf", dpi = 500, width = 2, height = 2.5)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "R_PF_ROI" = "#014421",
    "L_PF_ROI" = "#014421"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_motor.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('R_PF_ROI', 'L_PF_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#001a08") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_PF.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 2.5, height = 2)




###  6a  ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_6a_ROI', 'L_6a_ROI')) %>%
  group_by(lunaid, visitno, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#27763d"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#27763d")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 7, color = "#27763d"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_6a_ROI.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)


lm.model <- lmerTest::lmer("absBestError ~ ACW_50*age + epoch + (1|lunaid)", data = ACW_behav_avg)

summary(lm.model)

pred <- (ggpredict(lm.model, terms=c('ACW_50[all]', 'age[14, 20, 26]')))



lunaize(
  ggplot() +
    # Raw data points
    geom_point(
      data = ACW_behav_avg,
      aes(x = ACW_50, y = absBestError, shape=epoch),
      color = "#001a08",
      alpha = 0.25,
      size = 0.6
    ) + 
    geom_line(data = ACW_behav_avg,
              aes(x = ACW_50, y = absBestError, shape=epoch, group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
    # Model-predicted confidence intervals
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = conf.low, ymax = conf.high,
          fill = group),
      alpha = 0.15,
      color = NA
    ) +
    # Interaction lines
    geom_line(
      data = pred,
      aes(x = x, y = predicted, color = group),
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      name = "Age (years)"
    ) +
    scale_fill_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      guide = "none"
    )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank()
  ) +
  labs(
    x = "ACW (s)",
    y = "Accuracy (degrees)",
    color = "Age (years)"
  )




ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_6a_interaction.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)



jn <- interactions::johnson_neyman(
  lm.model,
  pred = ACW_50,
  modx = age
)

jn$plot +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  labs(
    title = NULL)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_6a_johnson.pdf", device = "pdf", dpi = 500, width = 1.7, height = 1.9)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "L_6a_ROI" = "#27763d"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_6a.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('R_6a_ROI', 'L_6a_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#27763d") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_6a.pdf", device = grDevices::cairo_pdf, dpi = 500,  width = 2.1, height = 1.4)





###  L_i6-8_ROI ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('L_i6-8_ROI', 'R_i6-8_ROI')) %>%
  group_by(lunaid, visitno, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#014421"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#014421")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 7, color = "#014421"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_i6-8_ROI.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)

lm.model <- lmerTest::lmer("absBestError ~ ACW_50*age + epoch + (1|lunaid)", data = ACW_behav_avg)

summary(lm.model)

pred <- (ggpredict(lm.model, terms=c('ACW_50[all]', 'age[14, 20, 26]')))



lunaize(
  ggplot() +
    # Raw data points
    geom_point(
      data = ACW_behav_avg,
      aes(x = ACW_50, y = absBestError, shape=epoch),
      color = "#001a08",
      alpha = 0.25,
      size = 0.6
    ) + 
    geom_line(data = ACW_behav_avg,
              aes(x = ACW_50, y = absBestError, shape=epoch, group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
    # Model-predicted confidence intervals
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = conf.low, ymax = conf.high,
          fill = group),
      alpha = 0.15,
      color = NA
    ) +
    # Interaction lines
    geom_line(
      data = pred,
      aes(x = x, y = predicted, color = group),
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      name = "Age (years)"
    ) +
    scale_fill_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      guide = "none"
    )) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      axis.line = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2),
      panel.grid = element_blank()
    ) +
    labs(
      x = "ACW (s)",
      y = "Accuracy (degrees)",
      color = "Age (years)"
    )




ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_LFL_interaction.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)


jn <- interactions::johnson_neyman(
  lm.model,
  pred = ACW_50,
  modx = age
)

jn$plot +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  labs(
    title = NULL)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_i6_8_johnson.pdf", device = "pdf", dpi = 500, width = 1.7, height = 1.9)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "L_i6-8_ROI" = "#014421"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_i6-8.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('L_i6-8_ROI', 'R_i6-8_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#014421") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_L_i6-8_ROI.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 2.1, height = 1.4)





###  L_IP2_ROI ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('L_IP2_ROI', 'R_IP2_ROI')) %>%
  group_by(lunaid, visitno, epoch) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#27763d"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#27763d")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 8, color = "#27763d"),
    axis.title.x = element_text(size = 8, color = c("black")),
    axis.title.y = element_text(size = 8, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_IP2_ROI.pdf", device = "pdf", dpi = 500, width = 2, height = 2.5)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "L_IP2_ROI" = "#27763d"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_IP2.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('L_IP2_ROI', 'L_IP2_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#001a08") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title.x = element_text(size = 7, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 7, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_IP2_ROI.pdf", device = grDevices::cairo_pdf, dpi = 500, width = 2.5, height = 2)





###  PEF ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('L_PEF_ROI', 'R_PEF_ROI')) %>%
  group_by(lunaid, visitno, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    absBestError = mean(absBestError, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(ACW_50 <0.05)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = absBestError, color = "#90b8a6"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#90b8a6")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 7, color = "#90b8a6"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Accuracy (degs)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_PEF_ROI.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)

lm.model <- lmerTest::lmer("absBestError ~ ACW_50*age + epoch + (1|lunaid)", data = ACW_behav_avg)

summary(lm.model)

pred <- (ggpredict(lm.model, terms=c('ACW_50[all]', 'age[14, 20, 26]')))



lunaize(
  ggplot() +
    # Raw data points
    geom_point(
      data = ACW_behav_avg,
      aes(x = ACW_50, y = absBestError, shape=epoch),
      color = "#001a08",
      alpha = 0.25,
      size = 0.6
    ) + 
    geom_line(data = ACW_behav_avg,
              aes(x = ACW_50, y = absBestError, shape=epoch, group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
    # Model-predicted confidence intervals
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = conf.low, ymax = conf.high,
          fill = group),
      alpha = 0.15,
      color = NA
    ) +
    # Interaction lines
    geom_line(
      data = pred,
      aes(x = x, y = predicted, color = group),
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      name = "Age (years)"
    ) +
    scale_fill_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      guide = "none"
    )) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank()
  ) +
  labs(
    x = "ACW (s)",
    y = "Accuracy (degrees)",
    color = "Age (years)"
  )




ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_PEF_interaction.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)



jn <- interactions::johnson_neyman(
  lm.model,
  pred = ACW_50,
  modx = age
)

jn$plot +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  labs(
    title = NULL)

ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_PEF_johnson.pdf", device = "pdf", dpi = 500, width = 1.7, height = 1.9)


ACW_behav_ROIs <- ACW_behav %>% merge(., glasser_rois, by = 'orig_parcelname', all=T) 


ggplot() +
  theme_void() +
  geom_brain(
    data = ACW_behav_ROIs,
    atlas = glasser,
    side = "lateral",
    hemi = "left",
    aes(fill = orig_parcelname)
  ) +
  scale_fill_manual(values = c(
    "L_PEF_ROI" = "#90b8a6"
  ), na.value = "white") +
  theme(
    legend.position = "none")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/brainPlots_PEF.pdf", device = "pdf", dpi = 500, width = 18, height = 5)

#### Neurosynth ----

neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>% dplyr::select(-regionName) %>% dplyr::select(-timing) #neurosynth term meta-analytic scores for 123 terms present in both NeuroSynth and the Cognitive Atlas, ordered in lh --> rh glasser order
neurosynth.termlist <- names(neurosynth.terms) #list of term names
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] #lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] #rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

sigRegions <- interaction_combined %>% filter(significant == TRUE) %>% filter(orig_parcelname %in% c('L_PEF_ROI', 'L_PEF_ROI'))

neurosynth.terms.a2 <- neurosynth.terms[neurosynth.terms$label %in% sigRegions$label,] %>% select(-label) #neurosynth z-scores in the  significant regions
neurosynth.terms.a2.weights <- colMeans(neurosynth.terms.a2) %>% as.data.frame() %>% set_names("term.mean")
neurosynth.terms.a2.weights$term <- row.names(neurosynth.terms.a2.weights)

neurosynth.terms.a2.weights %>% slice_max(n = 5, order_by = term.mean) %>%
  ggplot(., aes(x = term.mean, y = term)) +
  geom_segment(aes(y = reorder(term, term.mean), yend = term, x = 0, xend =
                     term.mean), color = "gray75", linewidth = .5) +
  geom_point(size = 1.5, color = "#90b8a6") +
  xlab("\nNeurosynth Term Score (z)") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6, family = "Arial", color = c("black")),
    axis.text.y = element_text(angle = 20, hjust = 1),
    axis.title.x = element_text(size = 6, family ="Arial", color = c("black")),
    axis.title.y = element_text(size = 6, family ="Arial", color = c("black")),
    axis.line = element_line(linewidth = .2),
    axis.ticks = element_line(linewidth = .2))



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/neuroSynth_PEF_ROI.pdf", device = grDevices::cairo_pdf, dpi = 500,  width = 2.1, height = 1.4)



# Absolute Best Error Variability ----

regions <- unique(ACW_behav$orig_parcelname)
ACW_behav$lunaid <- as.factor(ACW_behav$lunaid)
ACW_behav$epoch <- as.factor(ACW_behav$epoch)

# Initialize lists to store results
maineffect_list <- list()
interactions_list <- list()

# Loop through each combination of epoch and region
for (loc in regions) {
  
  lm.stats <- lm.model.function(input.df = ACW_behav, region = loc, dep_var = 'absBestError_sd', id_var = 'ACW_50')
  
  lm.stats.df <- data.frame(ACWestimate = lm.stats[2,1], ACWpValue = lm.stats[2,5], sexEstimate = lm.stats[4,1], sexpValue= lm.stats[4,5], epochEstimate = lm.stats[5,1], epochpValue = lm.stats[5,5]) %>% 
    mutate(orig_parcelname = loc)
  
  # Store results
  maineffect_list[[loc]] <- lm.stats.df

  
  
  lmInteraction.stats <- lm.interaction.model.function(input.df = ACW_behav, region = loc, dep_var = 'absBestError_sd', id_var = 'ACW_50')
  
  lmInteraction.stats.df <- data.frame(estimate = lmInteraction.stats[6,1], pValue = lmInteraction.stats[6,5])%>% 
    mutate(orig_parcelname = loc)
  
  
  interactions_list[[loc]] <- lmInteraction.stats.df
  
}


saveRDS(maineffect_list, file = "acw_behav_maineffects_linear_absBestError_sd.rds")
saveRDS(interactions_list, file = "acw_behav_interactioneffects_linear_absBestError_sd.rds")



maineffect_list <- readRDS("acw_behav_maineffects_linear_absBestError_sd.rds")
interactions_list <- readRDS("acw_behav_interactioneffects_linear_absBestError_sd.rds")


lmstats_combined <- bind_rows(maineffect_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

interaction_combined <- bind_rows(interactions_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )



glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

## run when remote ------------------------------------------------------------
# glasser <- glasser()
# glasser$data$hemi <- factor(glasser$data$hemi, levels = c('right','left'))
# glasser$data$side <- factor(glasser$data$side, levels = c('medial','lateral'))
# 

lmstats_combined <- lmstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)
interaction_combined <- interaction_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

## run when remote ------------------------------------------------------------
# lmstats_combined <- lmstats_combined %>%
#   mutate(label = paste0(label, "_ROI"))
# 
# interaction_combined <- interaction_combined %>%
#   mutate(label = paste0(label, "_ROI"))


dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#ffffff"   # white highlight
)

## main effect -----

lmstats_combined <- lmstats_combined %>%
  mutate(
    correctedP = p.adjust(ACWpValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate)
  ) +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )





## interaction -----

interaction_combined <- interaction_combined %>%
  mutate(
    correctedP = p.adjust(pValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate)
  ) +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )






# Latency ----
regions <- unique(ACW_behav$orig_parcelname)
ACW_behav$lunaid <- as.factor(ACW_behav$lunaid)
ACW_behav$epoch <- as.factor(ACW_behav$epoch)

# Initialize lists to store results
maineffect_list <- list()
interactions_list <- list()

# Loop through each combination of epoch and region
for (loc in regions) {
  
  lm.stats <- lm.model.function(input.df = ACW_behav, region = loc, dep_var = 'mgsLatency', id_var = 'ACW_50')
  
  lm.stats.df <- data.frame(ACWestimate = lm.stats[2,1], ACWpValue = lm.stats[2,5], sexEstimate = lm.stats[4,1], sexpValue= lm.stats[4,5], epochEstimate = lm.stats[5,1], epochpValue = lm.stats[5,5]) %>% 
    mutate(orig_parcelname = loc)
  
  # Store results
  maineffect_list[[loc]] <- lm.stats.df
  
  
  lmInteraction.stats <- lm.interaction.model.function(input.df = ACW_behav, region = loc, dep_var = 'mgsLatency', id_var = 'ACW_50')
  
  lmInteraction.stats.df <- data.frame(estimate = lmInteraction.stats[6,1], pValue = lmInteraction.stats[6,5])%>% 
    mutate(orig_parcelname = loc)
  
  
  interactions_list[[loc]] <- lmInteraction.stats.df
  
}


saveRDS(maineffect_list, file = "acw_behav_maineffects_linear_mgsLatency.rds")
saveRDS(interactions_list, file = "acw_behav_interactioneffects_linear_mgsLatency.rds")

maineffect_list <- readRDS("acw_behav_maineffects_linear_mgsLatency.rds")
interactions_list <- readRDS("acw_behav_interactioneffects_linear_mgsLatency.rds")


lmstats_combined <- bind_rows(maineffect_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

interaction_combined <- bind_rows(interactions_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )



glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

### run when remote ------------------------------------------------------------
# glasser <- glasser()
# glasser$data$hemi <- factor(glasser$data$hemi, levels = c('right','left'))
# glasser$data$side <- factor(glasser$data$side, levels = c('medial','lateral'))


lmstats_combined <- lmstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)
interaction_combined <- interaction_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

### run when remote ------------------------------------------------------------
# lmstats_combined <- lmstats_combined %>%
#   mutate(label = paste0(label, "_ROI"))
# 
# interaction_combined <- interaction_combined %>%
#   mutate(label = paste0(label, "_ROI"))


dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#ffffff"   # white highlight
)

## main effect -----

lmstats_combined <- lmstats_combined %>%
  mutate(
    correctedP = p.adjust(ACWpValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate)
  ) +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


## interaction ----

interaction_combined <- interaction_combined %>%
  mutate(
    correctedP = p.adjust(pValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate)
  ) +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) 






# Latency Var ----
regions <- unique(ACW_behav$orig_parcelname)
ACW_behav$lunaid <- as.factor(ACW_behav$lunaid)
ACW_behav$epoch <- as.factor(ACW_behav$epoch)

# Initialize lists to store results
maineffect_list <- list()
interactions_list <- list()

# Loop through each combination of epoch and region
for (loc in regions) {
  
  lm.stats <- lm.model.function(input.df = ACW_behav, region = loc, dep_var = 'mgsLatency_sd', id_var = 'ACW_50')
  
  lm.stats.df <- data.frame(ACWestimate = lm.stats[2,1], ACWpValue = lm.stats[2,5], sexEstimate = lm.stats[4,1], sexpValue= lm.stats[4,5], epochEstimate = lm.stats[5,1], epochpValue = lm.stats[5,5]) %>% 
    mutate(orig_parcelname = loc)
  
  # Store results
  maineffect_list[[loc]] <- lm.stats.df
  
  lmInteraction.stats <- lm.interaction.model.function(input.df = ACW_behav, region = loc, dep_var = 'mgsLatency_sd', id_var = 'ACW_50')
  
  lmInteraction.stats.df <- data.frame(estimate = lmInteraction.stats[6,1], pValue = lmInteraction.stats[6,5])%>% 
    mutate(orig_parcelname = loc)
  
  
  interactions_list[[loc]] <- lmInteraction.stats.df
  
}


saveRDS(maineffect_list, file = "acw_behav_maineffects_linear_mgsLatency_sd.rds")
saveRDS(interactions_list, file = "acw_behav_interactioneffects_linear_mgsLatency_sd.rds")

maineffect_list <- readRDS("acw_behav_maineffects_linear_mgsLatency_sd.rds")
interactions_list <- readRDS("acw_behav_interactioneffects_linear_mgsLatency_sd.rds")


lmstats_combined <- bind_rows(maineffect_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )

interaction_combined <- bind_rows(interactions_list, .id = "region") %>%
  separate(
    region,
    into = c("orig_parcelname"),
    sep = "_",
    extra = "merge",
    remove = TRUE
  )



glasser_rois <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/glasser360_regionlist.csv') %>%
  mutate(roi = gsub('_ROI', '', orig_parcelname)) %>%
  mutate(roi = gsub('L_', 'Left_', roi)) %>%
  mutate(roi = gsub('R_', 'Right_', roi))

## run when remote ------------------------------------------------------------
# glasser <- glasser()
# glasser$data$hemi <- factor(glasser$data$hemi, levels = c('right','left'))
# glasser$data$side <- factor(glasser$data$side, levels = c('medial','lateral'))


lmstats_combined <- lmstats_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)
interaction_combined <- interaction_combined %>% merge(., glasser_rois, by = 'orig_parcelname', all=T)

## run when remote ------------------------------------------------------------
# lmstats_combined <- lmstats_combined %>%
#   mutate(label = paste0(label, "_ROI"))
# 
# interaction_combined <- interaction_combined %>%
#   mutate(label = paste0(label, "_ROI"))


dramatic_palette <- c(
  "#001a08",  # near-black green
  "#012b0d",  # original dark forest
  "#014421",  # deep emerald
  "#27763d",  # original medium green
  "#4fa36a",  # brighter transition
  "#90b8a6",  # original sage
  "#d8efe2",  # pale mint
  "#ffffff"   # white highlight
)


## main effect ---- 

lmstats_combined <- lmstats_combined %>%
  mutate(
    correctedP = p.adjust(ACWpValue, method = "fdr"),
    significant = correctedP < 0.05)


ggplot() +
  theme_void() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate)
  ) +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = ACWestimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = lmstats_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, color = "#014421"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


## interaction ----

interaction_combined <- interaction_combined %>%
  mutate(
    correctedP = p.adjust(pValue, method = "fdr"),
    significant = correctedP < 0.05)

ggplot() +
  theme_void() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate)
  ) +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = estimate),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_gradientn(
    colours =  rev(dramatic_palette)) + 
  new_scale_fill() +
  geom_brain(
    data = interaction_combined,
    atlas = glasser,
    side = "lateral",
    aes(fill = significant),
    colour = "gray50",
    size = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = alpha("white", 0),
      "FALSE" = "gray80"
    ),
    na.value = alpha("white", 0)
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) 

## DLPFC ----

ACW_behav_avg <- ACW_behav %>%
  filter(orig_parcelname %in% c('R_a9-46v_ROI', 'L_a9-46v_ROI')) %>%
  group_by(lunaid, visitno, epoch, age) %>%
  summarise(
    ACW_50 = mean(ACW_50, na.rm = TRUE),
    mgsLatency_sd = mean(mgsLatency_sd, na.rm = TRUE),
    .groups = "drop"
  )%>% filter(ACW_50 <0.06)

lunaize(ggplot(data = ACW_behav_avg, aes(x = ACW_50, y = mgsLatency_sd, color = "#001a08"))) +
  geom_point(aes(shape = epoch), alpha = 0.3, size = 0.5) + 
  geom_line(aes(group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
  geom_smooth(aes(group=1), method='lm', alpha=0.3, linewidth = 0.8) + 
  scale_color_manual(values = c("#001a08")) + 
  scale_linetype_manual(values = c(
    "Rest" = "dashed",
    "Task" = "solid")) +
  (theme_minimal(base_family = "Arial")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = c("black")),
    plot.title = element_text(size = 7, color = "#001a08"),
    axis.title.x = element_text(size = 7, color = c("black")),
    axis.title.y = element_text(size = 7, color = c("black")),
    axis.line = element_line(linewidth = .2), 
    axis.ticks = element_line(linewidth = .2),
    strip.background = element_rect('white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  ggtitle('') +
  xlab("ACW (s)") +
  ylab("Lat Var (s)")


ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_dlpfc.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2)

lm.model <- lmerTest::lmer("mgsLatency_sd ~ ACW_50*age + epoch + (1|lunaid)", data = ACW_behav_avg)

summary(lm.model)

pred <- (ggpredict(lm.model, terms=c('ACW_50[all]', 'age[14, 20, 26]')))



lunaize(
  ggplot() +
    # Raw data points
    geom_point(
      data = ACW_behav_avg,
      aes(x = ACW_50, y = mgsLatency_sd, shape=epoch),
      color = "#001a08",
      alpha = 0.25,
      size = 0.6
    ) + 
    geom_line(data = ACW_behav_avg,
              aes(x = ACW_50, y = mgsLatency_sd, shape=epoch, group= interaction(lunaid, epoch)), alpha = 0.2, linewidth=0.1) +
    # Model-predicted confidence intervals
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = conf.low, ymax = conf.high,
          fill = group),
      alpha = 0.15,
      color = NA
    ) +
    # Interaction lines
    geom_line(
      data = pred,
      aes(x = x, y = predicted, color = group),
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      name = "Age (years)"
    ) +
    scale_fill_manual(
      values = c("#001a08", "#27763d", "#90b8a6"),
      guide = "none"
    )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank()
  ) +
  labs(
    x = "ACW (s)",
    y = "Accuracy (degrees)",
    color = "Age (years)"
  )



ggsave(filename = "/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/Figures/Figure4/acc_dlpfc_interaction.pdf", device = "pdf", dpi = 500, width = 1.5, height = 2.3)





