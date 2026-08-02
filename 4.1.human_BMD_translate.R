# this function runs the PBPK model get the deposition fraction after 24 hours 
# and then combine it together with the in vitro deposition fraction got from python 
# try to get the deposited amount before calculating the benchmark dose, instead, we calculate the 
# human equivalent dose first and then get the bmd calculation
library(stringr)
library(dplyr)
library(ggplot2)
source("with_stan/helper_function.R")

folder = "/Users/wuji/work/code/codo_v2/R/results/human/bmd_results/"


#----------1. PBPK model kinetic indicators ---------
#predicted_paras = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/pars_T_tot.csv")
#preditced_PK = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/mlr/PK_results.csv")
Body_weight = 73
F_W_Liver  = 0.0257  # Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401
preditced_human_PK = read.csv("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/PK_results_ratio.csv")
# load the processed effective concentration value
load("results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData")

#------2. read the in vitro DG dosimetry model results_vitro-------
vitro_dosimetry_df = read.csv("/Users/wuji/work/code/codo_v2/calculation_in_vitro/human/human_vitro_results_2025-10-14_21_33.csv")

vitro_dosimetry_df$Fraction_deposited_end 
# as for this file, we calculate for every concentration, so we directly use them all


# -------- 3. match the in vitro dosimetry model results with PROast results ---
name_mapping <- c(
  "SiO2 NPs" = "Silica",
  "SiO2" = "Silica",
  "Amorphous SiO2" = "Silica",
  "Graphene oxide" = "GO",
  "TiO2" = "TiO2",
  "Fe2O3" = "Iron Oxide"
)
vitro_dosimetry_df <- vitro_dosimetry_df %>%
  mutate(Substance = recode(`Substance`, !!!name_mapping))

vitro_dosimetry_df = shorten_cell_types(vitro_dosimetry_df,"cell_type")

combined_human_df_filtered<- combined_human_df_filtered %>%
  mutate(NP = recode(`NP`, !!!name_mapping))

matches <- lapply(seq_len(nrow(combined_human_df_filtered)), function(i) {
  row <- combined_human_df_filtered[i, ]
  subset(vitro_dosimetry_df,
         assay == row$assay &
           (Substance == row$NP_core | Substance == row$NP) &
           (Agglomerate_diameter_nm == row$HD) &
           In.vitro.media == row$media_type &
           Primary_particle_diameter_nm == row$Average_size_nm &
           cell_type == row$Cell_type &
           End_time_h == row$exposure_time)
})
length(matches) # check if it is one to one match

# Function to check the condition for each dataframe
check_unique <- function(lst) {
  sapply(lst, function(df) {
    col_name <- "Fraction_deposited_end"
    if (col_name %in% names(df)) {
      rounded_values <- round(unique(df[[col_name]]), digits = 8)  
      unique_vals = unique(rounded_values)
      if (length(unique_vals) == 1) {
        TRUE
      } else {
        unique_vals  # Return the unique values if not all the same
      }
    } else {
      FALSE  # Handle missing column as FALSE
      
    }
  })
}

# Apply the function to your list of dataframes
result <- check_unique(matches)

combined_human_df_filtered$Fraction_deposited_end <- sapply(matches, function(match_df) {
  if (nrow(match_df) > 0) {
    unique(match_df$Fraction_deposited_end)  # take first match
  } else {
    NA  # if no match, set to NA
  }
})
sum(is.na(combined_human_df_filtered$Fraction_deposited_end))

combined_human_df_filtered$Fraction_deposited_end = 1
combined_human_df_filtered$conc_deposited_end_median <- 
  combined_human_df_filtered$ec10_median * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q5  <- combined_human_df_filtered$ec10_q5  * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q25 <- combined_human_df_filtered$ec10_q25 * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q75 <- combined_human_df_filtered$ec10_q75 * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q95 <- combined_human_df_filtered$ec10_q95 * combined_human_df_filtered$Fraction_deposited_end


#--------------3. merge the in vitro model results_vitro with the PBPK model -------

library(dplyr)

combined_human_df_filtered <- combined_human_df_filtered %>%
  mutate(HD = ifelse(is.na(HD), Average_size_nm, HD))

preditced_human_PK<- preditced_human_PK %>%
  mutate(Hydrodynamic_Size = ifelse(is.na(Hydrodynamic_Size), Size, Hydrodynamic_Size))

# Loop through each "Substance name" in final_toxic_df
for (i in seq_len(nrow(combined_human_df_filtered))) {
  DE_L_id_g = 0
  
  # Extract the relevant data for the current row
  NP_core <- combined_human_df_filtered$`NP`[i]
  HD <- as.numeric(combined_human_df_filtered$`HD`[i])

  if (NP_core == "Graphene quantum dots (GQD)"){
    NP_core = "GO"
  }
  if (NP_core == "rGO"){
    NP_core = "GO"
  }
  # Retrieve the injected dose from final_toxic_df for the current row
  vitro_internal_conc_median <- combined_human_df_filtered$conc_deposited_end_median[i] # unit as mg/cm3
  vitro_internal_conc_q5     <- combined_human_df_filtered$conc_deposited_end_q5[i]       # unit: mg/cm³
  vitro_internal_conc_q25    <- combined_human_df_filtered$conc_deposited_end_q25[i]      # unit: mg/cm³
  vitro_internal_conc_q75    <- combined_human_df_filtered$conc_deposited_end_q75[i]      # unit: mg/cm³
  vitro_internal_conc_q95    <- combined_human_df_filtered$conc_deposited_end_q95[i]      # unit: mg/cm³
  
  density_liver = 1 #g/cm3
  
  
  # Filter AUC_result_all to only rows with matching NP_core 
  filtered_PK <- preditced_human_PK[grepl(NP_core, preditced_human_PK$NP_core), ]
  
  n_mc <- 10000
  set.seed(42)
  
  if (nrow(filtered_PK) > 0) {
    closest_row <- filtered_PK[which.min(abs(filtered_PK$Hydrodynamic_Size - HD)), ]
  } else {
    print(paste("No match found for NP_core:", NP_core, "instead find based on diameter"))
    closest_row <- preditced_human_PK[which.min(abs(preditced_human_PK$Hydrodynamic_Size - HD)), ]
  }
  
  # --- Extract DE quantiles ---
  DE_L_median  <- closest_row$human_DE_L_median
  DE_L_q5    <- closest_row$human_DE_L_q5
  DE_L_q95    <- closest_row$human_DE_L_q95
  
  DE24_L_median <- closest_row$human_DE24_L_median
  DE24_L_q5   <- closest_row$human_DE24_L_q5
  DE24_L_q95   <- closest_row$human_DE24_L_q95
  
  # --- Helper: Monte Carlo IVIVE ---
  mc_ivive <- function(bmd_median, bmd_q5, bmd_q95,
                       de_median, de_q5, de_q95,
                       liver_weight, density_liver, body_weight,
                       n_mc = 100000) {
    
    # --- BMD10 samples ---
    if (is.na(bmd_median) || bmd_median == 0) {
      return(list(quantiles = rep(NA, 5), 
                  var_bmd = NA, var_de = NA, 
                  pct_bmd = NA, pct_de = NA))
    }
    
    if (is.na(bmd_q5) || is.na(bmd_q95) || bmd_q5 == bmd_q95 || bmd_q5 == bmd_median) {
      bmd_samples <- rep(bmd_median, n_mc)
    } else {
      bmd_logsd <- (log(bmd_q95) - log(bmd_q5)) / (2 * 1.96)
      if (bmd_logsd <= 0) bmd_logsd <- 0.01
      bmd_samples <- rlnorm(n_mc, log(bmd_median), bmd_logsd)
    }
    
    var_bmd <- var(log10(bmd_samples))
    

    if (is.na(de_median) || de_median == 0) {
      return(list(quantiles = rep(NA, 5), 
                  var_bmd = var_bmd, var_de = NA, 
                  pct_bmd = NA, pct_de = NA))
    }
    
    # --- DE samples ---
    if (is.na(de_q5) || is.na(de_q95) || de_q5 == de_q95 || de_q5 == de_median) {
      de_samples <- rep(de_median, n_mc)
    } else {
      de_logsd <- (log(de_q95) - log(de_q5)) / (2 * 1.96)
      if (de_logsd <= 0) de_logsd <- 0.01
      de_samples <- rlnorm(n_mc, log(de_median), de_logsd)
    }
    
    var_de  <- var(log10(de_samples))
    pct_bmd <- var_bmd / (var_bmd + var_de) * 100
    pct_de  <- var_de / (var_bmd + var_de) * 100
    
    # IVIVE
    dose_samples <- (bmd_samples * liver_weight / density_liver) / 
      (de_samples * liver_weight / 100)
    dose_mg_kg   <- dose_samples / (body_weight * 1000)
    
    list(
      quantiles = quantile(dose_mg_kg, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)),
      var_bmd   = var_bmd,
      var_de    = var_de,
      pct_bmd   = pct_bmd,
      pct_de    = pct_de
    )
  }
  liver_weight <- F_W_Liver * Body_weight * 1000  # gram
  
  # --- DE at max time ---
  result <- mc_ivive(
    bmd_median = vitro_internal_conc_median,
    bmd_q5     = vitro_internal_conc_q5,
    bmd_q95    = vitro_internal_conc_q95,
    de_median  = DE_L_median,
    de_q5      = DE_L_q5,
    de_q95     = DE_L_q95,
    liver_weight = liver_weight,
    density_liver = density_liver,
    body_weight = Body_weight,
    n_mc = n_mc
  )
  
  dose_q <- result$quantiles
  
  # --- DE at 24h ---
  result_24 <- mc_ivive(
    bmd_median = vitro_internal_conc_median,
    bmd_q5     = vitro_internal_conc_q5,
    bmd_q95    = vitro_internal_conc_q95,
    de_median  = DE24_L_median,
    de_q5      = DE24_L_q5,
    de_q95     = DE24_L_q95,
    liver_weight = liver_weight,
    density_liver = density_liver,
    body_weight = Body_weight,
    n_mc = n_mc
  )
  
  dose_24_q <- result_24$quantiles
  
  # --- Save results ---
  if (!any(is.na(dose_q))) {
    combined_human_df_filtered[i, "Estimated_Dose_q5_mg_kg_BMD"]   <- dose_q[1]
    combined_human_df_filtered[i, "Estimated_Dose_q25_mg_kg_BMD"]    <- dose_q[2]
    combined_human_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- dose_q[3]
    combined_human_df_filtered[i, "Estimated_Dose_q75_mg_kg_BMD"]    <- dose_q[4]
    combined_human_df_filtered[i, "Estimated_Dose_q95_mg_kg_BMD"]   <- dose_q[5]
    
    combined_human_df_filtered[i, "Estimated_Dose_24_q5_mg_kg_BMD"]   <- dose_24_q[1]
    combined_human_df_filtered[i, "Estimated_Dose_24_q25_mg_kg_BMD"]    <- dose_24_q[2]
    combined_human_df_filtered[i, "Estimated_Dose_24_median_mg_kg_BMD"] <- dose_24_q[3]
    combined_human_df_filtered[i, "Estimated_Dose_24_q75_mg_kg_BMD"]    <- dose_24_q[4]
    combined_human_df_filtered[i, "Estimated_Dose_24_q95_mg_kg_BMD"]   <- dose_24_q[5]
    
    combined_human_df_filtered[i, "DE_L_median"]   <- DE_L_median
    combined_human_df_filtered[i, "DE24_L_median"] <- DE24_L_median
    combined_human_df_filtered[i, "closest_row_id"] <- closest_row$id
    
    combined_human_df_filtered[i, "var_log_bmd10"] <- result$var_bmd
    combined_human_df_filtered[i, "var_log_de"]    <- result$var_de
    combined_human_df_filtered[i, "pct_from_bmd10"] <- result$pct_bmd
    combined_human_df_filtered[i, "pct_from_de"]    <- result$pct_de
    
    combined_human_df_filtered[i, "var_log_bmd10_24"]   <- result_24$var_bmd
    combined_human_df_filtered[i, "var_log_de_24"]      <- result_24$var_de
    combined_human_df_filtered[i, "pct_from_bmd10_24"]  <- result_24$pct_bmd
    combined_human_df_filtered[i, "pct_from_de_24"]     <- result_24$pct_de
    
  } else {
    combined_human_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- NA
    combined_human_df_filtered[i, "DE_L_median"]   <- DE_L_median
    combined_human_df_filtered[i, "DE24_L_median"] <- DE24_L_median
    combined_human_df_filtered[i, "closest_row_id"] <- closest_row$id
  }
  
}

summary(combined_human_df_filtered$var_log_bmd10)
summary(combined_human_df_filtered$var_log_bmd10_24)

summary(combined_human_df_filtered$pct_from_bmd10)
summary(combined_human_df_filtered$pct_from_de)
summary(combined_human_df_filtered$pct_from_bmd10_24)
summary(combined_human_df_filtered$pct_from_de_24)



library(dplyr)

table_variance_decomposition <- combined_human_df_filtered %>%
  group_by(NP) %>%
  summarise(
    n = n(),
    bmd_median = round(median(pct_from_bmd10, na.rm = TRUE), 1),
    bmd_min    = round(min(pct_from_bmd10, na.rm = TRUE), 1),
    bmd_max    = round(max(pct_from_bmd10, na.rm = TRUE), 1),
    de_median  = round(median(pct_from_de, na.rm = TRUE), 1),
    de_min     = round(min(pct_from_de, na.rm = TRUE), 1),
    de_max     = round(max(pct_from_de, na.rm = TRUE), 1)
  ) %>%
  mutate(
    bmd_col = paste0(bmd_median, " [", bmd_min, "–", bmd_max, "]"),
    de_col  = paste0(de_median, " [", de_min, "–", de_max, "]")
  )


all_row <- combined_human_df_filtered %>%
  summarise(
    NP = "All",
    n = n(),
    bmd_col = paste0(
      round(median(pct_from_bmd10, na.rm = TRUE), 1), " [",
      round(min(pct_from_bmd10, na.rm = TRUE), 1), "–",
      round(max(pct_from_bmd10, na.rm = TRUE), 1), "]"
    ),
    de_col = paste0(
      round(median(pct_from_de, na.rm = TRUE), 1), " [",
      round(min(pct_from_de, na.rm = TRUE), 1), "–",
      round(max(pct_from_de, na.rm = TRUE), 1), "]"
    )
  )

table_variance_decomposition <- bind_rows(table_variance_decomposition, all_row)

write.csv(table_variance_decomposition, 
          "results_undissolved/human/variance_decomposition_summary_human.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
#-------------- 4. plot the results-----
plot_folder = "results_undissolved/human/plots/"
results_vitro = combined_human_df_filtered
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/"
save(results_vitro, file = file.path(stanoutdir,"translated_dose_human_all.RData"))

#load(file.path(stanoutdir,"translated_dose_human_all.RData"))
results_vitro_insoluble <- results_vitro[!(results_vitro$Substance %in% c("Ag", "ZnO")), ]
results_vitro_insoluble <- results_vitro[(results_vitro$NP %in% c("rGO","GO", "Iron Oxide","Silica","TiO2","Au")), ]
write.csv(results_vitro_insoluble, file.path(stanoutdir,"translated_dose_human_all.csv"), row.names = FALSE)

# split into different assay category
results_vitro$assay_category <- ifelse(
  grepl("viability|Cytotoxicity", results_vitro$`assay`, ignore.case = TRUE), 
  "viability", 
  ifelse(
    grepl("ROS|TNF|IL|DCFDA", results_vitro$`assay`, ignore.case = TRUE), 
    "oxidative_stress", 
    results_vitro$`assay`
  )
)
# split into different assay category
results_vitro_insoluble$assay_category <- ifelse(
  grepl("viability|Cytotoxicity", results_vitro_insoluble$`assay`, ignore.case = TRUE), 
  "viability", 
  ifelse(
    grepl("ROS|TNF|IL|DCFDA", results_vitro_insoluble$`assay`, ignore.case = TRUE), 
    "oxidative_stress", 
    results_vitro_insoluble$`assay`
  )
)



library(ggplot2)
results_vitro_insoluble$general_experiment_type[
  results_vitro_insoluble$general_experiment_type == "Biochemical"
] <- "Metabolic Cell Stress"


vitro_toxicity_plot = ggplot(
  results_vitro_insoluble,
  aes(x = NP, y = Estimated_Dose_median_mg_kg_BMD,fill=general_experiment_type)
) +
  # 5–95% interval (thin)
  geom_linerange(
    aes(ymin = Estimated_Dose_q5_mg_kg_BMD,
        ymax = Estimated_Dose_q95_mg_kg_BMD),
        linewidth = 1, color = "#808080"
  ) +
  # 25–75% interval (thicker)
  geom_linerange(
    aes(ymin = Estimated_Dose_q25_mg_kg_BMD,
        ymax = Estimated_Dose_q75_mg_kg_BMD),
    linewidth = 3, color = "#808080",alpha= 0.7
  ) +
  # Median point
  geom_point(
    shape = 21, size = 4, stroke = 1.5
  ) +
  scale_y_log10(limits = c(1e-2, 1e3)) +
  scale_fill_manual(values = c( "#BF5A2D","#8BAFCD", "#6C3687","#DAB039"),
                    guide  = guide_legend(order = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(hjust = 1, size = 12, face = "bold",color = "black"),
    axis.text.y  = element_text(size = 12, face = "bold",color = "black"),
    axis.title.x = element_text(size = 14, face = "bold",color = "black"),
    axis.title.y = element_text(size = 14, face = "bold",color = "black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1.5)
  ) +
  labs(
    #title = "In vivo Injection Dosage Leading to Liver Cell Concentrations Reaching EC10",
    x = "Nanoparticle",
    y = bquote(bold("Translated")~ bolditalic("in vivo") ~ bold("BMD"[10]^IVIVE~"(mg/kg BW)"))
  ) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21), 
                        byrow = TRUE, order = 1, title = "Experiment type")
  )
vitro_toxicity_plot
#ggsave(paste0(plot_folder,"translated_dose/insoluble/estimated_dose_vitro_insoluble.pdf"),
#       plot = vitro_toxicity_plot, width = 8, height = 5)


range_df <- results_vitro_insoluble %>%
  group_by(NP) %>%
  summarise(
    Estimated_Dose_median_mg_kg_BMD = median(
      Estimated_Dose_median_mg_kg_BMD, na.rm = TRUE
    ),
    
    # overall 5–95% range for each NP
    Estimated_Dose_q5_mg_kg_BMD = min(
      Estimated_Dose_q5_mg_kg_BMD, na.rm = TRUE
    ),
    Estimated_Dose_q95_mg_kg_BMD = max(
      Estimated_Dose_q95_mg_kg_BMD, na.rm = TRUE
    ),
    
    # overall 25–75% range for each NP
    Estimated_Dose_q25_mg_kg_BMD = min(
      Estimated_Dose_q25_mg_kg_BMD, na.rm = TRUE
    ),
    Estimated_Dose_q75_mg_kg_BMD = max(
      Estimated_Dose_q75_mg_kg_BMD, na.rm = TRUE
    ),
    
    .groups = "drop"
  )

plot_human_bmd_ivive_combined = ggplot() +
  geom_boxplot(
    data = range_df,
    aes(
      x = NP,
      ymin = Estimated_Dose_q5_mg_kg_BMD,
      lower = Estimated_Dose_q25_mg_kg_BMD,
      middle = Estimated_Dose_median_mg_kg_BMD,
      upper = Estimated_Dose_q75_mg_kg_BMD,
      ymax = Estimated_Dose_q95_mg_kg_BMD
    ),
    stat = "identity",
    width = 0.35,
    linewidth = 0.8,
    color = "#808080",
    fill = "#F0EBE1"

  ) +
  
  # original points, colored by experiment type
  geom_point(
    data = results_vitro_insoluble,
    aes(
      x = NP,
      y = Estimated_Dose_median_mg_kg_BMD,
      fill = general_experiment_type
    ),
    shape = 21,
    size = 4,
    stroke = 1.5,
    color = "black",
    position = position_jitter(width = 0.12, height = 0)
  ) +
  
  scale_y_log10() +
  scale_fill_manual(
    values = c("#BF5A2D", "#8BAFCD", "#6C3687", "#DAB039"),
    guide = guide_legend(order = 1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      hjust = 1, size = 12, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 12, face = "bold", color = "black"
    ),
    axis.title.x = element_text(
      size = 14, face = "bold", color = "black"
    ),
    axis.title.y = element_text(
      size = 14, face = "bold", color = "black"
    ),
    panel.grid.major = element_line(
      color = "gray80", linetype = "dashed"
    ),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(
      fill = "white", color = "black", linewidth = 1.5
    )
  ) +
  labs(
    x = "Nanoparticle",
    y = bquote(
      bold("Translated") ~ bolditalic("in vivo") ~
        bold("BMD"[10]^IVIVE ~ "(mg/kg BW)")
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21),
      byrow = TRUE,
      order = 1,
      title = "Experiment type"
    )
  )
plot_human_bmd_ivive_combined
ggsave(paste0(plot_folder,"translated_dose/insoluble/estimated_dose_vitro_insoluble_combined.pdf"),
       plot = plot_human_bmd_ivive_combined, width = 8, height = 5)

library(dplyr)
library(ggplot2)


point_df <- results_vitro_insoluble %>%
  group_by(NP) %>%
  mutate(
    x_label = paste0(.data[["experiment_type"]], " No.", .data[["study"]])
  ) %>%
  ungroup()

plot_human_bmd_ivive_separate = ggplot(
  point_df,
  aes(x = x_label)
) +
  geom_boxplot(
    aes(
      ymin   = Estimated_Dose_q5_mg_kg_BMD,
      lower  = Estimated_Dose_q25_mg_kg_BMD,
      middle = Estimated_Dose_median_mg_kg_BMD,
      upper  = Estimated_Dose_q75_mg_kg_BMD,
      ymax   = Estimated_Dose_q95_mg_kg_BMD,
      fill   = general_experiment_type
    ),
    stat = "identity",
    width = 0.45,
    linewidth = 0.3,
    color = "gray30"
  ) +
  
  facet_wrap(~ NP, scales = "free_x") +
  scale_y_log10() +
  scale_fill_manual(
    values = c("#BF5A2D", "#8BAFCD", "#6C3687", "#DAB039"),
    guide = guide_legend(order = 1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x = element_text(
      angle = 90, hjust = 1, size = 9, face = "bold", color = "black"
    ),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(
      fill = "white", color = "black", linewidth = 1
    )
  ) +
  labs(
    x = "Experiment",
    y = bquote(
      bold("Translated") ~ bolditalic("in vivo") ~
        bold("BMD"[10]^IVIVE ~ "(mg/kg BW)")
    ),
    fill = "Experiment type"
  )
plot_human_bmd_ivive_separate
ggsave(paste0(plot_folder,"translated_dose/insoluble/estimated_dose_vitro_insoluble_separate.pdf"),
       plot = plot_human_bmd_ivive_separate, width = 15, height = 10)

# ------- plot with human clinical results comparison -----

results_agg <- results_vitro_insoluble %>%
  group_by(NP) %>%
  summarise(
    ymin   = min(Estimated_Dose_q5_mg_kg_BMD, na.rm = TRUE),
    lower  = min(Estimated_Dose_q25_mg_kg_BMD, na.rm = TRUE),
    middle = median(Estimated_Dose_median_mg_kg_BMD, na.rm = TRUE),
    upper  = max(Estimated_Dose_q75_mg_kg_BMD, na.rm = TRUE),
    ymax   = max(Estimated_Dose_q95_mg_kg_BMD, na.rm = TRUE)
  ) %>%
  filter(is.finite(middle))
auroshell_row <- data.frame(
  NP     = "AuroShell",
  ymin   = 22,
  lower  = 22,
  middle = (22 + 36) / 2,
  upper  = 36,
  ymax   = 36
)
results_agg <- bind_rows(results_agg, auroshell_row)
results_agg$NP <- factor(results_agg$NP,
                          levels = c("AuroShell", 
                                     setdiff(unique(results_agg$NP), "AuroShell")))

agg_vitro_toxicity_plot = ggplot(results_agg, aes(x = NP)) +
  geom_boxplot(
    aes(lower = lower, upper = upper, middle = middle,
        ymin = ymin, ymax = ymax,
        fill = ifelse(NP == "AuroShell", "Clinical reference", "IVIVE-translated dose")),
    stat = "identity",
    color = "black",
    width = 0.5, alpha = 0.7
  ) +
  scale_y_log10() +
  scale_fill_manual(values = c("Clinical reference" = "#DAB039",
                               "IVIVE-translated dose" = "#8BAFCD")) +
  #coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text  = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1.5)
  ) +
  labs(
    y = bquote(bold("Translated") ~ bolditalic("in vivo") ~ bold("BMD"[10]^IVIVE ~ "(mg/kg BW)")),
    x = "Nanoparticle"
  )
  
  
  
agg_vitro_toxicity_plot
ggsave(paste0(plot_folder,"translated_dose/insoluble/estimated_dose_vitro_insoluble_agg.pdf"),
       plot = agg_vitro_toxicity_plot, width = 8, height = 5)

