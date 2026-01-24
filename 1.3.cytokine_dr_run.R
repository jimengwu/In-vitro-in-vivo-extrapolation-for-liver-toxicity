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
outputdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/plots/bmd_human_cytokine_stan"
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
      title = paste0("Study ", i, ": cytokine results"),
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
  ggsave(
    file.path(outputdir, sprintf("study_%d_curve.pdf", i)),
    plot = p, width = 6, height = 5,device = cairo_pdf) 
}


#============== calculate using not only the dose respnse curve but the NOEC value ==============

# Initialize result dataframe
min_signif_results <- data.frame(
  study = numeric(),
  min_signif_conc = character(),  # use character to allow "no_toxic"
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
  threshold <- 0.90  # adjust based on biological relevance
  signif_idx <- which(diffs > threshold)
  if (length(signif_idx) > 0) {
    min_signif_conc <- min(obs_df$x[signif_idx])
  } else {
    min_signif_conc <- "no_toxic"
  }
  cat("Minimal significantly different concentration:", min_signif_conc, "\n")
  # Add to results dataframe
  min_signif_results <- rbind(
    min_signif_results,
    data.frame(
      study = i,
      min_signif_conc = as.character(min_signif_conc),
      stringsAsFactors = FALSE
    )
  )
}

# for mouse
save(min_signif_results, file = file.path(stanoutdir,"cytokine_mouse_ec5_draws_all.RData"))
# for human
save(min_signif_results, file = file.path(stanoutdir,"cytokine_human_ec5_draws_all.RData"))

