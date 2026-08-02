# ==============================================================================
#
# TITLE:   Threshold derivation for cytokine / inflammatory experiments
#
# PURPOSE: Takes the cleaned, per-experiment in vitro data and derives the
#          toxicity threshold for the cytokine endpoint class (TNF-a, IL-6,
#          IL-8, IL-1b).
#          These readouts are frequently non-monotonic (e.g. bell-shaped) or
#          threshold-triggered and are usually reported as fold-change versus
#          control, so no sigmoidal Hill model is imposed. Instead, a LOEC-based
#          surrogate is used: the threshold is the lowest tested concentration
#          producing a change from the negative control beyond the predefined
#          benchmark (90% change in this study). This value substitutes for the
#          BMD10 in all downstream IVIVE steps.
#
# INPUT:   <cleaned per-experiment data set from 01_read_clean_invitro_data.R>
# OUTPUT:  <per-experiment threshold table (LOEC-based BMD10 surrogate) + plots>
#
# NOTE:    No curve fitting and no credible intervals here — thresholds are
#          restricted to the tested dose grid, which limits their resolution.
#
# AUTHOR:  Jimeng Wu          CREATED: 2026-04          
# ==============================================================================

library(cmdstanr)
library(ggplot2)
library(dplyr)
library(stringr)

# mouse cytokine dataset 
load("/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/ls_sub_case_mouse.RData")
outputdir = "results_undissolved/mouse/plots/bmd_mouse_cytokine_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/stan_results/bmd_mouse_cytokine_stan"

# human cytokine dataset
load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData")
outputdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/plots/bmd_human_cytokine_stan_new"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results/bmd_human_cytokine_stan"

#===================== start of the code =========================
ls_bmd <- ls_sub_case[sapply(ls_sub_case, nrow) > 1]
ls_bmd_cytokine <- ls_bmd[grep("cytokine", names(ls_bmd), 
                           ignore.case = TRUE)]

unique_tested_assays <- unique(unlist(lapply(ls_bmd_cytokine, function(x) x$`Tested assay`)))
unique_tested_assays


results <- data.frame(study = integer(), a = numeric(), b = numeric(),c = numeric(), 
                      d = numeric(), sigma = numeric(), ec5_median = numeric(), 
                      ec5_lower = numeric(), ec5_upper = numeric(),ec10_median = numeric(), 
                      ec10_lower = numeric(), ec10_upper = numeric())
ec5_draws_all <- list()
ec10_draws_all <- list()
param_draws_all <- list()
r_squared_df <- list()
for (i in seq_along(ls_bmd_cytokine)) {
  data_i <- ls_bmd_cytokine[[i]]
  p = ggplot() +
    geom_point(
      data = data_i,
      aes(x = `In vitro concentration`, y = Results),
      size = 2
    )+
    labs(
      #title = paste0("Case ", i, ": cytokine results"),
      title = stringr::str_wrap(paste0("Cytokine", " case ", i, ": ", 
                                       data_i$`Substance name`[1],"_", 
                                       data_i$`Cell type`[1]), width = 50),
      x = "In vitro concentration (μg/ml)",
      y =  data_i$`Tested assay` %>%
        unique() %>%
        str_remove("^Cytokine:\\s*") %>%
        str_replace("Interleukin-8 assay \\(IL-8\\)", "IL-8"))  +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 14,color = "black"),   # x-axis text size
      axis.title = element_text(size = 16, face = "bold",color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  p
  #ggsave(
  #  file.path(outputdir, sprintf("study_%d_curve.pdf", i)),
  #  plot = p, width = 6, height = 5,device = cairo_pdf) 
}


#============== calculate using not only the dose respnse curve but the NOEC value ==============

# Initialize result dataframe
min_signif_results <- data.frame(
  study = numeric(),
  lower_conc = character(),
  lower_fold = character(),
  min_signif_conc = character(),
  min_signif_conc_fold = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(ls_bmd_cytokine)) {
  print(i)
  data_i <- ls_bmd_cytokine[[i]]
  x_obs <- data_i$`In vitro concentration`
  y_obs <- data_i$Results
  obs_df <- data.frame(x = x_obs, y = y_obs)
  
  assay_type <- unique(data_i$`Tested assay`)
  print(assay_type)
  

  
  # Check if "fold of control" is in the assay type, if not, skip this experiment, needs preprocess for it
  if (!any(grepl("fold of control", assay_type, ignore.case = TRUE))) {
    obs_df$y <- (obs_df$y / min(obs_df$y, na.rm = TRUE))
  }else if (any(grepl("% of control", assay_type, ignore.case = TRUE))) {
    obs_df$y = obs_df$y/100
  }
    
  # Ensure a reference point (x = 0, y = 1) exists in obs_df
  if (!0 %in% obs_df$x) {
    obs_df <- rbind(obs_df, data.frame(x = 0, y = 1))
  }

  obs_df <- obs_df[order(obs_df$x), ]  # Sort by concentration (ascending)
  
  # extract the minimal concentration when it is significantly different with benchmark value which is y = 1
  
  # Simple threshold-based alternative if only 1 observation per concentration
  diffs <- obs_df$y - 1
  # --- Threshold: 2-fold (100% change) ---
  threshold <- 1
  signif_idx <- which(diffs > threshold)
  
  if (length(signif_idx) > 0) {
    first_pos <- min(signif_idx)
    loec <- obs_df$x[first_pos]
    loec_fold <- obs_df$y[first_pos]
    if (first_pos > 1) {
      lower_conc <- obs_df$x[first_pos - 1]
      lower_fold <- obs_df$y[first_pos - 1]
    } else {
      lower_conc <- 0
      lower_fold <- 1
    }
  } else {
    loec <- "no_toxic"
    loec_fold <- NA
    lower_conc <- NA
    lower_fold <- NA
  }
  
  # --- Flag whether LOEC changed ---
  min_signif_results <- rbind(
    min_signif_results,
    data.frame(
      study = i,
      lower_conc= as.character(lower_conc),
      lower_fold = as.character(lower_fold),
      min_signif_conc = as.character(loec),
      min_signif_conc_fold = as.character(loec_fold),
      stringsAsFactors = FALSE
    )
  )
}
# for mouse
write.csv(min_signif_results, file.path(stanoutdir,"cytokine_mouse_loec_draws_all.csv"), row.names = FALSE)
save(min_signif_results, file = file.path(stanoutdir,"cytokine_mouse_ec5_draws_all.RData"))

# for human
write.csv(min_signif_results, file.path(stanoutdir,"cytokine_human_loec_draws_all.csv"), row.names = FALSE)
save(min_signif_results, file = file.path(stanoutdir,"cytokine_human_ec5_draws_all.RData"))