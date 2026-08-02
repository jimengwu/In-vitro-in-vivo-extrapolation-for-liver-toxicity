# TITLE:   Subgroup analysis of BMD10 by experimental covariates
#          
#
# PURPOSE: Addresses the concern that the compiled in vitro dataset spans
#          multiple cell lines (HepG2, L-02, Kupffer, primary hepatocytes),
#          exposure durations (6-72 h) and assay types, which are otherwise
#          treated as interchangeable within each nanoparticle-assay category -
#          risking confusion between true biological sensitivity and
#          experimental artefacts.
#          The script stratifies the experiment-level BMD10 estimates by cell
#          line, exposure duration and assay type, and tests whether these
#          factors systematically shift the BMD10 distribution. Comparisons are
#          made within nanoparticle type where the sample size allows, and
#          group sizes are reported alongside each result so that
#          underpowered strata are not over-interpreted.
#
# INPUT:   <experiment-level BMD10 table with covariate annotation
#           (cell line, exposure time, assay, nanoparticle, species)>
# OUTPUT:  <stratified BMD10 summary tables>
#          <subgroup comparison figures (intended for the supporting information)>
#
# NOTE:    Descriptive / exploratory analysis. Several strata contain only a few
#          experiments, so the outcome is reported as an absence of a systematic
#          trend rather than as evidence of no effect.
#
# AUTHOR:  Jimeng Wu          CREATED: 2026-04          
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# 0. PREREQUISITES — adjust paths to match your project
# ==============================================================================
# Assumes combined_human_df_filtered and combined_mouse_df_filtered 
# are already loaded in your workspace (same objects used in your main analysis).
# If not, source your data-loading script here:
# source("your_data_loading_script.R")
load("results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData")
load("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_filtered_results.RData")

# Filter to non-dissolvable NPs (same filter as your main analysis)
human_df <- combined_human_df_filtered %>%
  filter(NP %in% c("Au", "TiO2", "SiO2", "GO", "rGO"))

mouse_df <- combined_mouse_df_filtered %>%
  filter(NP %in% c("Au", "TiO2", "SiO2", "GO", "rGO"))


# ==============================================================================
# 1. DATA OVERVIEW — how many experiments per cell line, exposure time, NP
# ==============================================================================

# --- Human ---
human_summary <- human_df %>%
  mutate(
    exposure_bin = case_when(
      exposure_time <= 24 ~ "≤24 h",
      exposure_time > 24 & exposure_time <= 48 ~ "25–48 h",
      exposure_time > 48 ~ ">48 h",
      TRUE ~ "Unknown"
    ),
    Average_size_nm = as.numeric(Average_size_nm),
    HD = as.numeric(HD)
  ) %>%
  group_by(NP, general_experiment_type, Cell_type, exposure_bin) %>%
  summarise(
    n_experiments = n(),
    size_median = round(median(Average_size_nm, na.rm = TRUE), 1),
    size_min    = round(min(Average_size_nm, na.rm = TRUE), 1),
    size_max    = round(max(Average_size_nm, na.rm = TRUE), 1),
    HD_median   = round(median(HD, na.rm = TRUE), 1),
    HD_min      = round(min(HD, na.rm = TRUE), 1),
    HD_max      = round(max(HD, na.rm = TRUE), 1),
    n_HD_available = sum(!is.na(HD)),
    BMD10_median = round(median(ec10_median, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(NP, general_experiment_type, Cell_type)
cat("\n=== Human: Experiment counts by NP / Assay / Cell type / Exposure bin ===\n")
print(human_summary)
write.csv(human_summary, "results_undissolved/sub_group_analysis/Table_S_human_BMD10_summary.csv", row.names = FALSE)

# --- Mouse ---
mouse_summary <- mouse_df %>%
  mutate(exposure_bin = case_when(
    exposure_time <= 24 ~ "≤24 h",
    exposure_time > 24 & exposure_time <= 48 ~ "25–48 h",
    exposure_time > 48 ~ ">48 h",
    TRUE ~ "Unknown"
  )) %>%
  count(NP, general_experiment_type, Cell_type, exposure_bin, name = "n_experiments") %>%
  arrange(NP, general_experiment_type, Cell_type)

cat("\n=== Mouse: Experiment counts by NP / Assay / Cell type / Exposure bin ===\n")
print(mouse_summary)

# ==============================================================================
# 2. STRATIFIED ANALYSIS — by Cell Line
# ==============================================================================

# --- Figure S_CellLine: BMD10 stratified by cell line, faceted by NP ---

# Human
p_cellline_human <- human_df %>%
  filter(!is.na(ec10_median), ec10_median > 0) %>%
  ggplot(aes(x = Cell_type, y = ec10_median, color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", ncol = 2) +
  labs(
    title = "Human in vitro BMD10 stratified by cell line",
    x = "Cell type",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("Figure_S_cellline_human.png", p_cellline_human,
       width = 10, height = 8, dpi = 300)
ggsave("Figure_S_cellline_human.pdf", p_cellline_human,
       width = 10, height = 8)

# Mouse
p_cellline_mouse <- mouse_df %>%
  filter(!is.na(ec10_median), ec10_median > 0) %>%
  ggplot(aes(x = Cell_type, y = ec10_median, color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", ncol = 2) +
  labs(
    title = "Mouse in vitro BMD10 stratified by cell line",
    x = "Cell type",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("Figure_S_cellline_mouse.png", p_cellline_mouse,
       width = 10, height = 8, dpi = 300)
ggsave("Figure_S_cellline_mouse.pdf", p_cellline_mouse,
       width = 10, height = 8)

# ==============================================================================
# 3. STRATIFIED ANALYSIS — by Exposure Duration
# ==============================================================================

# Create exposure duration bins
add_exposure_bin <- function(df) {
  df %>%
    mutate(
      exposure_time_num = as.numeric(exposure_time),
      exposure_bin = case_when(
        exposure_time_num <= 24 ~ "≤24 h",
        exposure_time_num > 24 & exposure_time_num <= 48 ~ "25–48 h",
        exposure_time_num > 48 ~ ">48 h",
        TRUE ~ "Unknown"
      ),
      exposure_bin = factor(exposure_bin, levels = c("≤24 h", "25–48 h", ">48 h", "Unknown"))
    )
}

human_df_exp <- add_exposure_bin(human_df)
mouse_df_exp <- add_exposure_bin(mouse_df)

# --- Figure S_ExposureTime: BMD10 stratified by exposure duration ---

# Human
p_exposure_human <- human_df_exp %>%
  filter(!is.na(ec10_median), ec10_median > 0, exposure_bin != "Unknown") %>%
  ggplot(aes(x = exposure_bin, y = ec10_median, color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", ncol = 2) +
  labs(
    title = "Human in vitro BMD10 stratified by exposure duration",
    x = "Exposure duration",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("Figure_S_exposure_human.png", p_exposure_human,
       width = 10, height = 8, dpi = 300)
ggsave("Figure_S_exposure_human.pdf", p_exposure_human,
       width = 10, height = 8)

# Mouse
p_exposure_mouse <- mouse_df_exp %>%
  filter(!is.na(ec10_median), ec10_median > 0, exposure_bin != "Unknown") %>%
  ggplot(aes(x = exposure_bin, y = ec10_median, color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", ncol = 2) +
  labs(
    title = "Mouse in vitro BMD10 stratified by exposure duration",
    x = "Exposure duration",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("Figure_S_exposure_mouse.png", p_exposure_mouse,
       width = 10, height = 8, dpi = 300)
ggsave("Figure_S_exposure_mouse.pdf", p_exposure_mouse,
       width = 10, height = 8)

# ==============================================================================
# 4. COMBINED PANEL FIGURE — Cell line + Exposure in one figure (for SI)
# ==============================================================================

# This creates a cleaner combined figure for the human dataset
# Panel A: by cell line, Panel B: by exposure duration

human_plot_df <- human_df_exp %>%
  filter(!is.na(ec10_median), ec10_median > 0)

# Panel A
pA <- human_plot_df %>%
  ggplot(aes(x = Cell_type, y = ec10_median, 
             color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", nrow = 1) +
  labs(
    subtitle = "(A) Stratified by cell line",
    x = "Cell type",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  )

# Panel B
pB <- human_plot_df %>%
  filter(exposure_bin != "Unknown") %>%
  ggplot(aes(x = exposure_bin, y = ec10_median, 
             color = general_experiment_type)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", nrow = 1) +
  labs(
    subtitle = "(B) Stratified by exposure duration",
    x = "Exposure duration",
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Combine using patchwork (install if needed: install.packages("patchwork"))
library(patchwork)
p_combined <- pA / pB + plot_layout(heights = c(1, 1.15))


ggsave("results_undissolved/sub_group_analysis/Figure_S_subgroup_combined_human.pdf", p_combined,
       width = 10, height = 6)

# ==============================================================================
# 5. SUMMARY STATISTICS TABLE — median BMD10 by subgroup
# ==============================================================================

# --- Table S: Summary statistics by cell line ---
table_cellline <- human_df %>%
  filter(!is.na(ec10_median), ec10_median > 0) %>%
  group_by(NP, general_experiment_type, Cell_type) %>%
  summarise(
    n = n(),
    median_BMD10 = round(median(ec10_median, na.rm = TRUE), 2),
    Q25 = round(quantile(ec10_median, 0.25, na.rm = TRUE), 2),
    Q75 = round(quantile(ec10_median, 0.75, na.rm = TRUE), 2),
    min_BMD10 = round(min(ec10_median, na.rm = TRUE), 2),
    max_BMD10 = round(max(ec10_median, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(NP, general_experiment_type, Cell_type)

cat("\n=== Table S: Human BMD10 summary by cell line ===\n")
print(table_cellline, n = Inf)
write.csv(table_cellline, "Table_S_BMD10_by_cellline.csv", row.names = FALSE)

# --- Table S: Summary statistics by exposure duration ---
table_exposure <- human_df_exp %>%
  filter(!is.na(ec10_median), ec10_median > 0, exposure_bin != "Unknown") %>%
  group_by(NP, general_experiment_type, exposure_bin) %>%
  summarise(
    n = n(),
    median_BMD10 = round(median(ec10_median, na.rm = TRUE), 2),
    Q25 = round(quantile(ec10_median, 0.25, na.rm = TRUE), 2),
    Q75 = round(quantile(ec10_median, 0.75, na.rm = TRUE), 2),
    min_BMD10 = round(min(ec10_median, na.rm = TRUE), 2),
    max_BMD10 = round(max(ec10_median, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(NP, general_experiment_type, exposure_bin)

cat("\n=== Table S: Human BMD10 summary by exposure duration ===\n")
print(table_exposure, n = Inf)
write.csv(table_exposure, "Table_S_BMD10_by_exposure.csv", row.names = FALSE)

# ==============================================================================
# 6. META-REGRESSION — which factors explain BMD10 variance?
# ==============================================================================

# Prepare the regression dataset
reg_df <- human_df_exp %>%
  filter(!is.na(ec10_median), ec10_median > 0) %>%
  mutate(
    log10_BMD10 = log10(ec10_median),
    log10_size  = log10(as.numeric(Average_size_nm)),
    exposure_time_num = as.numeric(exposure_time)
  ) %>%
  filter(!is.na(log10_BMD10))

# --- Model 1: Full model with all available covariates ---
# NOTE: adjust column names if your data uses different names
model_full <- lm(
  log10_BMD10 ~ NP + general_experiment_type + Cell_type + 
    exposure_time_num + log10_size,
  data = reg_df
)

cat("\n=== Meta-regression: Full model ===\n")
print(summary(model_full))

# --- Model 2: Material properties only (NP type + size) ---
model_material <- lm(
  log10_BMD10 ~ NP + log10_size,
  data = reg_df %>% filter(!is.na(log10_size))
)

cat("\n=== Meta-regression: Material properties only ===\n")
print(summary(model_material))

# --- Model 3: Experimental factors only (cell line + exposure time) ---
model_experimental <- lm(
  log10_BMD10 ~ Cell_type + exposure_time_num,
  data = reg_df
)

cat("\n=== Meta-regression: Experimental factors only ===\n")
print(summary(model_experimental))

# --- Compare R² values to see which factors explain more variance ---
cat("\n=== R² comparison ===\n")
cat(sprintf("Full model R²:                %.3f (adj: %.3f)\n",
            summary(model_full)$r.squared,
            summary(model_full)$adj.r.squared))
cat(sprintf("Material properties only R²:  %.3f (adj: %.3f)\n",
            summary(model_material)$r.squared,
            summary(model_material)$adj.r.squared))
cat(sprintf("Experimental factors only R²: %.3f (adj: %.3f)\n",
            summary(model_experimental)$r.squared,
            summary(model_experimental)$adj.r.squared))

# --- ANOVA decomposition of the full model ---
cat("\n=== ANOVA: Sequential sum of squares (Type I) ===\n")
print(anova(model_full))

# --- Export regression results as a table ---
coef_table <- broom::tidy(model_full, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

write.csv(coef_table, "results_undissolved/sub_group_analysis/Table_S_metaregression_coefficients.csv", row.names = FALSE)
cat("\n=== Regression coefficients saved to Table_S_metaregression_coefficients.csv ===\n")

# ==============================================================================
# 7. FOCUSED COMPARISON — HepG2 vs other cell types (Comment 5)
# ==============================================================================

# For the reviewer's specific concern about HepG2 bias
human_df_cellcomp <- human_df %>%
  filter(!is.na(ec10_median), ec10_median > 0) %>%
  mutate(cell_group = ifelse(Cell_type == "HepG2", "HepG2", "Other liver cells"))

# Summary by HepG2 vs others, within each NP
table_hepg2 <- human_df_cellcomp %>%
  group_by(NP, general_experiment_type, cell_group) %>%
  summarise(
    n = n(),
    median_BMD10 = round(median(ec10_median, na.rm = TRUE), 2),
    geometric_mean = round(10^mean(log10(ec10_median), na.rm = TRUE), 2),
    Q25 = round(quantile(ec10_median, 0.25, na.rm = TRUE), 2),
    Q75 = round(quantile(ec10_median, 0.75, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(NP, general_experiment_type, cell_group)

cat("\n=== HepG2 vs other liver cells — BMD10 comparison ===\n")
print(table_hepg2, n = Inf)
write.csv(table_hepg2, "Table_S_HepG2_vs_others.csv", row.names = FALSE)

# Plot: HepG2 vs others
p_hepg2 <- human_df_cellcomp %>%
  ggplot(aes(x = cell_group, y = ec10_median, color = general_experiment_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ NP, scales = "free_x", nrow = 1) +
  labs(
    title = "BMD10 comparison: HepG2 vs other liver cell models",
    x = NULL,
    y = expression(BMD[10] ~ "(µg/mL)"),
    color = "Assay category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("Figure_S_HepG2_comparison.png", p_hepg2,
       width = 12, height = 5, dpi = 300)
ggsave("Figure_S_HepG2_comparison.pdf", p_hepg2,
       width = 12, height = 5)

# ==============================================================================
# 8. OPTIONAL — Wilcoxon test for HepG2 vs others (per NP, where n allows)
# ==============================================================================

# Only run where both groups have ≥3 observations
test_results <- human_df_cellcomp %>%
  group_by(NP) %>%
  filter(n_distinct(cell_group) == 2) %>%
  summarise(
    n_HepG2 = sum(cell_group == "HepG2"),
    n_Other = sum(cell_group == "Other liver cells"),
    .groups = "drop"
  ) %>%
  filter(n_HepG2 >= 3 & n_Other >= 3)

if (nrow(test_results) > 0) {
  cat("\n=== Wilcoxon rank-sum tests: HepG2 vs Other (log10 BMD10) ===\n")
  for (np in test_results$NP) {
    sub <- human_df_cellcomp %>%
      filter(NP == np, !is.na(ec10_median), ec10_median > 0)
    
    wt <- wilcox.test(
      log10(ec10_median) ~ cell_group,
      data = sub,
      exact = FALSE
    )
    
    cat(sprintf(
      "%s: W = %.1f, p = %.4f (n_HepG2 = %d, n_Other = %d)\n",
      np, wt$statistic, wt$p.value,
      sum(sub$cell_group == "HepG2"),
      sum(sub$cell_group == "Other liver cells")
    ))
  }
} else {
  cat("\nInsufficient data for formal HepG2 vs Other statistical tests ",
      "(need ≥3 per group within a NP type).\n")
}

cat("\n=== Analysis complete. Outputs saved. ===\n")