

# --- Setup (same as before, done once) ---
neurosynth.terms <- read.csv("/Volumes/Hera/Projects/corticalmyelin_development/Maps/Neurosynth_terms/atl-glasser360_neurosynth.csv") %>%
  dplyr::select(-regionName) %>%
  dplyr::select(-timing) 
neurosynth.termlist <- names(neurosynth.terms)
neurosynth.terms$label[1:180] <- glasser_rois$label[181:360] # lh first
neurosynth.terms$label[181:360] <- glasser_rois$label[1:180] # rh
neurosynth.terms[,1:123] <- scale(neurosynth.terms[,1:123])

# --- All significant regions (no parcelname filter now) ---
sigRegions <- interaction_combined %>% 
  filter(significant == TRUE)

sig_labels <- unique(sigRegions$label)

# --- Loop through each significant region ---
region_top5_plots <- list()
region_top5_terms <- list()

for (lbl in sig_labels) {
  
  region_scores <- neurosynth.terms %>%
    filter(label == lbl) %>%
    dplyr::select(-label)
  
  if (nrow(region_scores) == 0) next  # skip if no match
  
  region_scores_long <- region_scores %>%
    tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "z_score")
  
  top5 <- region_scores_long %>%
    slice_max(n = 5, order_by = z_score)
  
  region_top5_terms[[lbl]] <- top5
  
  p <- ggplot(top5, aes(x = z_score, y = reorder(term, z_score))) +
    geom_segment(aes(y = reorder(term, z_score), yend = term, x = 0, xend = z_score),
                 color = "gray75", linewidth = .5) +
    geom_point(size = 1.5, color = "#001a08") +
    xlab("\nNeurosynth Term Score (z)") +
    ylab("") +
    ggtitle(lbl) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 8, family = "Arial", color = "black"),
      axis.text.y = element_text(angle = 0, hjust = 1),
      axis.title.x = element_text(size = 7, family = "Arial", color = "black"),
      axis.title.y = element_text(size = 7, family = "Arial", color = "black"),
      axis.line = element_line(linewidth = .2),
      axis.ticks = element_line(linewidth = .2)
    )
  
  region_top5_plots[[lbl]] <- p
}

# --- View a specific region's plot, e.g. ---
region_top5_plots[["lh_L_1"]]

region_top5_df <- dplyr::bind_rows(region_top5_terms, .id = "label") %>%
  dplyr::left_join(sigRegions %>% dplyr::select(label, orig_parcelname) %>% dplyr::distinct(), 
                   by = "label")
