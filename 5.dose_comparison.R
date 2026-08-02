#------------- read the mouse in vivo data -------
library(readxl)
library(dplyr)


in_vivo_df = read_excel("/Users/wuji/work/code/codo_v2/data/mouse_invivo.xlsx")
results_vivo = in_vivo_df[2:nrow(in_vivo_df),c("doi","material", "Particle primary diameter nm",
                                               "Type",
                                               "endpoint_type","endpoint_value",
                                               "injected_dose_acute(mg/kg)",
                                               "injected_dose_LD50(mg/kg)",
                                               "experiment type","Animal")]

results_vivo = unique(results_vivo)

results_vivo <- results_vivo %>%
  mutate(ec10_eq = case_when(
    endpoint_type %in% c("LOEC", "LOEC-") ~ endpoint_value / 2,
    endpoint_type == "LD50"               ~ endpoint_value / 10,
    TRUE                                  ~ endpoint_value
  ))

results_vivo_unique = results_vivo %>%
  distinct(endpoint_type, material, `Particle primary diameter nm`, doi,`Type`,`endpoint_type`,
           `endpoint_value`,`experiment type`, `ec10_eq`)
results_vivo_unique = results_vivo_unique %>%
  filter(material %in% c("GO", "Silica","TiO2"))


results_vivo_unique <- results_vivo_unique %>%
  mutate(
    `experiment type` = if_else(
      `experiment type` %in% c("mortality"),
      "mortality",
      "no_mortality"
    )
  )


stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/"
load(file = file.path(stanoutdir,"translated_dose_mouse_all.RData"))

combined_mouse_df_filtered <- combined_mouse_df_filtered[(combined_mouse_df_filtered$NP %in% c("GO", "Iron Oxide","Silica","TiO2","Au")), ]

combined_mouse_df_filtered

combined_df <- bind_rows(
  combined_mouse_df_filtered %>% 
    select(study,Paper_ID, NP, Estimated_Dose_median_mg_kg_BMD, ec10_median,
           general_experiment_type,Cell_type,Average_size_nm,HD) %>%
    rename(Dose = Estimated_Dose_median_mg_kg_BMD, Substance = NP) %>%
    mutate(
      Source = paste("in vitro", general_experiment_type, sep = "_"),
      Dose = as.numeric(Dose)
    ),
  
  results_vivo_unique %>%
    select(material, `ec10_eq`,`Particle primary diameter nm`,`endpoint_type`,
           `experiment type`) %>%
    rename(Dose = `ec10_eq`, Substance = material) %>%
    mutate(Source = "in vivo",
          #Source = `experiment type`,
           Dose = as.numeric(Dose))
)

# ------------ quantify the invitro in vivo alignment -----------------
# ---  fold difference ---
quant_comparison <- combined_df %>%
  mutate(source_group = ifelse(grepl("in vivo", Source), "vivo", "vitro")) %>%
  group_by(Substance, source_group) %>%
  summarise(
    n = n(),
    median_dose = median(Dose, na.rm = TRUE),
    gm = exp(mean(log(Dose), na.rm = TRUE)),  # geometric mean
    q25 = quantile(Dose, 0.25, na.rm = TRUE),
    q75 = quantile(Dose, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = source_group,
    values_from = c(n, median_dose, gm, q25, q75)
  ) %>%
  mutate(
    fold_diff_median = median_dose_vitro / median_dose_vivo,
    fold_diff_gm = gm_vitro / gm_vivo,
    log10_fold = log10(fold_diff_gm)
  )

print(quant_comparison)

quant_comparison %>% select(Substance, fold_diff_gm)

#  in vitro within in vivo range 10 folds
within_factor_summary <- overall %>%
  group_by(Substance) %>%
  summarise(
    vivo_min = min(Dose[source_group == "vivo"], na.rm = TRUE),
    vivo_max = max(Dose[source_group == "vivo"], na.rm = TRUE),
    n_vitro = sum(source_group == "vitro"),
    n_within_10x = sum(
      source_group == "vitro" &
        Dose >= vivo_min / 10 &
        Dose <= vivo_max * 10
    ),
    pct_within_10x = round(n_within_10x / n_vitro * 100, 1),
    .groups = "drop"
  )

print(within_factor_summary)

within_factor_summary %>% select(Substance, pct_within_10x)

# --- Overall log10 RMSD ---
overall_rmsd <- quant_comparison %>%
  filter(!is.na(gm_vitro) & !is.na(gm_vivo)) %>%
  summarise(
    log10_RMSD = sqrt(mean((log10(gm_vitro) - log10(gm_vivo))^2)),
    fold_deviation = 10^sqrt(mean((log10(gm_vitro) - log10(gm_vivo))^2))
  )

print(overall_rmsd)



# ------------ quantify the invitro in vivo alignment -----------------


plot_folder = "results_undissolved/mouse/plots/"
#plot_mouse_comparison(data=combined_df, 
#                      output_path= paste0(plot_folder, 
#                                          "translated_dose/insoluble/comparison_estimated_dose_mouse_insoluble_BMD2.pdf"),
#                      plot_title = "In Vivo Dosages vs. In Vitro-Translated BMD10")

library(ggplot2)

p <- ggplot(combined_df, aes(x = Source, y = Dose, fill = Source)) +
  geom_violin(
    color = "black", alpha = 0.6, trim = FALSE, width = 0.8
  ) +
  geom_jitter(
    width = 0.15, size = 1.5, alpha = 1, shape = 21, color = "black"
  ) +
  facet_wrap(~ Substance, nrow = 1, scales = "free_x") +
  scale_y_log10() +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(
    name   = "Data source",
    values = c(
      "in vitro_Cytotoxicity" = "#8BAFCD",
      "in vitro_Biochemical"  = "#DAB039",
      "in vitro_Cytokine"     = "#BF5A2D",
      "in vivo"               = "#d5d1d1"
    ),breaks = c("in vitro_Cytotoxicity", "in vitro_Biochemical", "in vitro_Cytokine", "in vivo"),
    labels = c(
      expression(italic("In vitro") ~ cytotoxicity),
      expression(italic("In vitro") ~ "metabolic cell stress"),
      expression(italic("In vitro") ~ cytokine),
      expression(italic("In vivo") ~ equivalent)
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(size = 16, face = "bold", color = "black"),
    strip.text         = element_text(size = 14, face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    legend.title       = element_text(color = "black", size = 14, face = "bold"),
    legend.text        = element_text(color = "black", size = 12, face = "bold"),
    panel.background   = element_rect(fill = "white", color = "black", linewidth = 1.5)
  ) +
  labs(y = bquote(bold(bolditalic("In vivo") ~ "BMD"[10] ~ "(mg/kg BW)")))
p
ggsave(paste0(plot_folder, 
              "translated_dose/insoluble/comparison_estimated_dose_mouse_insoluble_BMD_classified.pdf"),
       plot = p, width = 14, height = 4)

# ======================= plot with size comparison =============================
library(dplyr)
library(readr)

plot_df <- combined_df %>%
  mutate(
    diameter_nm = coalesce(
      parse_number(as.character(Average_size_nm)),
      parse_number(as.character(`Particle primary diameter nm`))
    ),
    diameter_class = case_when(
      is.na(diameter_nm) ~ "Unknown",
      diameter_nm > 50 ~ "> 50 nm",
      TRUE ~ "<= 50 nm"
    ),
    diameter_class = factor(
      diameter_class,
      levels = c("<= 50 nm", "> 50 nm", "Unknown")
    )
  )

p_size <- ggplot(plot_df, aes(x = Substance, y = Dose, fill = Source)) +
  
  geom_violin(
    position = position_dodge(width = 0.8),
    color = "black",
    alpha = 0.6,
    trim = FALSE,
    width = 0.8
  ) +
  
  geom_jitter(
    aes(
      shape = diameter_class,
      group = Source
    ),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.8
    ),
    size = 2.5,
    alpha = 1,
    color = "black"
  ) +
  
  scale_y_log10() +
  
  scale_fill_manual(
    name = "Data source",
    breaks = c(
      "in vitro_Cytotoxicity",
      "in vitro_Biochemical",
      "in vitro_Cytokine",
      "in vivo"
    ),
    values = c(
      "in vitro_Cytotoxicity" = "#8BAFCD",
      "in vitro_Biochemical"  = "#DAB039",
      "in vitro_Cytokine"     = "#BF5A2D",
      "in vivo"               = "#d5d1d1"
    ),
    labels = c(
      expression(italic("In vitro") ~ cytotoxicity),
      expression(italic("In vitro") ~ "metabolic stress"),
      expression(italic("In vitro") ~ cytokine),
      expression(italic("In vivo") ~ equivalent)
    )
  ) +
  
  scale_shape_manual(
    name = "Primary diameter",
    values = c(
      "<= 50 nm" = 4,
      "> 50 nm"  = 21,
      "Unknown"  = 23
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    
    legend.title = element_text(color = "black", size = 14, face = "bold"),
    legend.text  = element_text(color = "black", size = 12, face = "bold"),
    
    panel.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 1.5
    )
  ) +
  
  labs(
    x = "Nanoparticle core",
    y = bquote(bold(bolditalic("In vivo") ~ "BMD"[10] ~ "(mg/kg BW)"))
  ) +guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "black",
        size = 3,
        stroke = 1.1,
        alpha = 1
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        color = "black",
        size = 3,
        stroke = 1.1
      )
    )
  )

p_size
ggsave(paste0(plot_folder, 
              "translated_dose/insoluble/comparison_estimated_dose_mouse_insoluble_BMD_classified_size.pdf"),
       plot = p_size, width = 10, height = 5)
#- ------ VIOLIN PLOT for THE IN VITRO VS IN VIVO IN MOUSE

# --- build spacers (k = number of empty slots between substances) ---
combined_df <- combined_df %>%
  mutate(Group = ifelse(grepl("in vitro", Source), "in vitro", "in vivo"))

k <- 2
subs <- if (is.factor(combined_df$Substance)) levels(combined_df$Substance) else unique(combined_df$Substance)

lvls <- unlist(lapply(seq_along(subs), function(i) {
  s <- subs[i]
  if (i < length(subs)) c(s, paste0(s, "_sp", seq_len(k))) else s  # no spacers after last
}))
combined_df$Substance_spaced <- factor(combined_df$Substance, levels = lvls)

# optional split lines AFTER each block (not after the last since no spacer_k there)
xsplit <- which(grepl(paste0("_sp", k, "$"), levels(combined_df$Substance_spaced))) + 0.5

pd <- position_dodge(width = 5)

p_violin = ggplot(combined_df, aes(x = Substance_spaced, y = Dose, fill = Group)) +
  #geom_vline(xintercept = xsplit, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  #geom_violin(trim = FALSE, color = "black", width = 3, 
  #            alpha = 0.6, position = pd) +
  geom_boxplot(aes(group = interaction(Substance_spaced, Group)),
               width = 1.5, outlier.shape = NA, alpha = 0.6, color = "black",
               position = position_dodge(width = 2)) +
  geom_point(aes(fill = Group), color = "black", stroke = 0.3, shape = 21,
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 2),
             size = 2, show.legend = FALSE) +
  scale_y_log10(limit=c(1e-3,1e4)) +
  # show only real Substance labels; no trailing space on the right border
  scale_x_discrete(
    drop = FALSE,
    breaks = subs, labels = subs,
    expand = expansion(mult = c(0, 0), add = c(0.03, 0))  # left margin only; right = 0
  ) +
  #scale_fill_manual(values = c("in vitro" = "#F6DFD6", "in vivo" = "#B6B3D6")) +
  scale_fill_manual(
    name   = "Data source",
    values = c("in vitro" = "#ECB66C", "in vivo" = "#BFDFD2"),
    labels = c(expression(italic("in vitro")), expression(italic("in vivo")))
  )+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold",color="black"),
    axis.text.y  = element_text(size = 12, face = "bold",color="black"),
    axis.title.x = element_text(size = 14, face = "bold",color="black"),
    axis.title.y = element_text(size = 14, face = "bold",color="black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray80", linetype = "dashed"),
    legend.title     = element_text(color="black", size=12, face="bold"),
    legend.text      = element_text(color="black", size=12,face="bold"),
    panel.background = element_rect(fill = "white", color = "black", size = 1.2)
  ) +
  labs(
    x = "Nanoparticle",
    y = "Dose (mg/kg)"
  )
p_violin
ggsave(paste0(plot_folder, 
              "translated_dose/insoluble/comparison_estimated_dose_mouse_insoluble_BMD.pdf"), 
       plot = p_violin, width = 6, height = 4)


# ---------------- for human ---------

stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/"
load(file = file.path(stanoutdir,"translated_dose_human_all.RData"))

results_vitro <- results_vitro[(results_vitro$NP 
                                                          %in% c("GO", "Iron Oxide","Silica",
                                                                 "TiO2","Au")), ]
plot_human_dose = results_vitro %>% 
  select(NP, Estimated_Dose_median_mg_kg_BMD, general_experiment_type,Average_size_nm) %>%
  rename(Dose = Estimated_Dose_median_mg_kg_BMD, Substance = NP) %>%
  mutate(
    Source = paste("in vitro", general_experiment_type, sep = "_"),
    Dose = as.numeric(Dose)
  )

ggplot(plot_human_dose, aes(x = Substance, y = Dose)) +
  #geom_vline(xintercept = xsplit, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  #geom_violin(trim = FALSE, color = "black", width = 3, 
  #            alpha = 0.6, position = pd) +
  geom_boxplot(width = 1.5, outlier.shape = NA, alpha = 0.6, color = "black",
               position = position_dodge(width = 2)) +

  scale_y_log10(limit=c(1e-5,1e5)) +
  # show only real Substance labels; no trailing space on the right border

  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold",color="black"),
    axis.text.y  = element_text(size = 12, face = "bold",color="black"),
    axis.title.x = element_text(size = 14, face = "bold",color="black"),
    axis.title.y = element_text(size = 14, face = "bold",color="black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray80", linetype = "dashed"),
    legend.title     = element_text(color="black", size=12, face="bold"),
    legend.text      = element_text(color="black", size=12,face="bold"),
    panel.background = element_rect(fill = "white", color = "black", size = 1.2)
  ) +
  labs(
    x = "Nanoparticle",
    y = "Dose (mg/kg)"
  )

