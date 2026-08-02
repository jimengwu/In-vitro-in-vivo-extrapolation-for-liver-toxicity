# this function runs the PBPK model get the deposition fraction after 24 hours 
# and then combine it together with the in vitro deposition fraction got from python 
# try to get the deposited amount before calculating the benchmark dose, instead, we calculate the 
# human equivalent dose first and then get the BMD calculation


library(dplyr)
library(ggplot2)
library(stringr)

folder = "/Users/wuji/work/code/codo_v2/R/results/mouse/bmd_results/"

#-------------0. setting-------
Body_weight = 0.02  # kg, Bodyweight (Davies 1993)
F_W_Liver  = 0.055  # (Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401)
F_W_Brain  = 0.017 # Unitless, Fractional brain tissue (Brown 1997)
F_W_Lung   = 0.007  # Unitless, Fractional lung tissue (Brown 1997)
F_W_Kidney = 0.017  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Spleen = 0.005  # Unitless®, Fractional spleen tissue (Brown 1997)


#----------1. PBPK model kinetic indicators ---------
predicted_paras = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/pars_T_tot.csv")
preditced_PK = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/mlr/PK_results.csv")
merged_PK <- merge( preditced_PK,predicted_paras,by.x = "id", by.y = "Folder", all = TRUE)

mouse_AUC = read.csv("~/work/code/Mouse-general-PBPK/plots/auc_results.csv")
merged_PK <- merge( mouse_AUC,preditced_PK ,by.x = "id", by.y = "id", all = TRUE)


#------2. read the in vitro DG dosimetry model results_vitro-------
vitro_dosimetry_df = read.csv("/Users/wuji/work/code/codo_v2/calculation_in_vitro/mouse/mouse_vitro_results_2025-10-14_14_52.csv")
vitro_dosimetry_df$Fraction_deposited_end  # as for this file, we calculate for every concentration, so we directly use them all
vitro_dosimetry_df <- vitro_dosimetry_df %>%
  mutate(cell_type = str_replace(cell_type, "^Hepa1-6$", "Hepa 1-6"))

final_toxic_df = vitro_dosimetry_df

# ------- 2.1 load BMD results calculated from in vitro data directly -----

load("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_filtered_results.RData")


# -------- 3. match the in vitro dosimetry model results with PROast results ---
# Define a mapping from old names to new names
name_mapping <- c(
  "SiO2 NPs" = "Silica",
  "SiO2" = "Silica",
  "Amorphous SiO2" = "Silica",
  "Graphene oxide" = "GO",
  "TiO2" = "TiO2",
  "Fe2O3" = "Iron Oxide"
  #"ZnO" = "Zinc Oxide",
  #"Ag" = "Silver"
)
vitro_dosimetry_df <- vitro_dosimetry_df %>%
  mutate(Substance = recode(`Substance`, !!!name_mapping))


combined_mouse_df_filtered<- combined_mouse_df_filtered %>%
  mutate(NP = recode(`NP`, !!!name_mapping))
# Match each row from combined_mouse_df_filtered to vitro_dosimetry_df

matches_mouse <- lapply(seq_len(nrow(combined_mouse_df_filtered)), function(i) {
  row <- combined_mouse_df_filtered[i, ]
  subset(vitro_dosimetry_df,
         assay == row$assay &
           (Substance == row$NP_core | Substance == row$NP) &
           Agglomerate_diameter_nm == row$HD &
           In.vitro.media == row$media_type &
           Primary_particle_diameter_nm == row$Average_size_nm &
           cell_type == row$Cell_type &
           End_time_h == row$exposure_time)
})



# Extract the first matching 'Fraction_deposited_end' value (or NA)
fraction_deposited <- sapply(matches_mouse, function(df) {
  if (nrow(df) > 0) unique(df$Fraction_deposited_end) else NA_real_
})
fraction_deposited = 1
# Add to the original data frame
combined_mouse_df_filtered$Fraction_deposited_end <- fraction_deposited

# Calculate deposited concentration
combined_mouse_df_filtered$conc_deposited_end_median <- 
  combined_mouse_df_filtered$ec10_median * combined_mouse_df_filtered$Fraction_deposited_end
combined_mouse_df_filtered$conc_deposited_end_q5  <- combined_mouse_df_filtered$ec10_q5  * combined_mouse_df_filtered$Fraction_deposited_end
combined_mouse_df_filtered$conc_deposited_end_q25 <- combined_mouse_df_filtered$ec10_q25 * combined_mouse_df_filtered$Fraction_deposited_end
combined_mouse_df_filtered$conc_deposited_end_q75 <- combined_mouse_df_filtered$ec10_q75 * combined_mouse_df_filtered$Fraction_deposited_end
combined_mouse_df_filtered$conc_deposited_end_q95 <- combined_mouse_df_filtered$ec10_q95 * combined_mouse_df_filtered$Fraction_deposited_end

#--------------3. merge the in vitro model results_vitro with the PBPK model -------

# Loop through each "Substance name" in final_toxic_df
for (i in seq_len(nrow(combined_mouse_df_filtered))) {
  DE_L_id_g = 0
  
  # Extract the relevant data for the current row
  NP_core <- combined_mouse_df_filtered$`NP`[i]
  HD <- as.numeric(combined_mouse_df_filtered$`HD`[i])

  
  # Retrieve the injected dose from final_toxic_df for the current row
  vitro_internal_conc_median <- combined_mouse_df_filtered$conc_deposited_end_median[i] # unit as mg/cm3
  vitro_internal_conc_q5     <- combined_mouse_df_filtered$conc_deposited_end_q5[i]       # unit: mg/cm³
  vitro_internal_conc_q25    <- combined_mouse_df_filtered$conc_deposited_end_q25[i]      # unit: mg/cm³
  vitro_internal_conc_q75    <- combined_mouse_df_filtered$conc_deposited_end_q75[i]      # unit: mg/cm³
  vitro_internal_conc_q95    <- combined_mouse_df_filtered$conc_deposited_end_q95[i]      # unit: mg/cm³
  
  density_liver = 1 #g/cm3
  
  
  # Filter merged_PK to only rows with matching NP_core #todo: ZnO is not found,
  filtered_PK <- merged_PK[grepl(NP_core, merged_PK$NP_core), ]
  
  if (nrow(filtered_PK) > 0) {
    # Find the row with the closest Hydrodynamic_Size
    closest_row <- filtered_PK[which.min(abs(filtered_PK$Hydrodynamic_Size - HD)), ]
    # Extract the corresponding DE_L_id_g
    DE_L_id_g <- closest_row$DE_L_id_g
    DE24_L_id_g <- closest_row$DE24_L_id_g
    
  } else {
    print(paste("No match found for NP_core:", NP_core,"instead find based on diameter"))
    closest_row <- merged_PK[which.min(abs(merged_PK$Hydrodynamic_Size - HD)), ]
    # Extract the corresponding DE_L_id_g
    DE_L_id_g <- closest_row$DE_L_id_g
    DE24_L_id_g <- closest_row$DE24_L_id_g
  }
  liver_weight = F_W_Liver * Body_weight * 1000 # gram as unit
  
  Dose_q5     <- (vitro_internal_conc_q5     * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q25    <- (vitro_internal_conc_q25    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_median <- (vitro_internal_conc_median * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q75    <- (vitro_internal_conc_q75    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q95    <- (vitro_internal_conc_q95    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  
  if (length(Dose_median) > 0) {
    combined_mouse_df_filtered[i, "Estimated_Dose_q5_mg_kg_BMD"]     <- Dose_q5     / (Body_weight * 1000)  # µg → mg/kg
    combined_mouse_df_filtered[i, "Estimated_Dose_q25_mg_kg_BMD"]    <- Dose_q25    / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- Dose_median / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_q75_mg_kg_BMD"]    <- Dose_q75    / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_q95_mg_kg_BMD"]    <- Dose_q95    / (Body_weight * 1000)
    combined_mouse_df_filtered[i,"closest_row_id"] <- closest_row$id
    combined_mouse_df_filtered[i, "DE_L_id_g"] <- DE_L_id_g
    combined_mouse_df_filtered[i, "DE24_L_id_g"] <- DE24_L_id_g
    } else {
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- NA
    combined_mouse_df_filtered[i,"closest_row_id"] <- closest_row$id
    combined_mouse_df_filtered[i, "DE_L_id_g"] <- DE_L_id_g
    combined_mouse_df_filtered[i, "DE24_L_id_g"] <- DE24_L_id_g
  }
  
}


for (i in seq_len(nrow(combined_mouse_df_filtered))) {
  DE_L_id_g = 0
  
  NP_core <- combined_mouse_df_filtered$`NP`[i]
  HD <- as.numeric(combined_mouse_df_filtered$`HD`[i])
  
  vitro_internal_conc_median <- combined_mouse_df_filtered$conc_deposited_end_median[i]
  vitro_internal_conc_q5     <- combined_mouse_df_filtered$conc_deposited_end_q5[i]
  vitro_internal_conc_q25    <- combined_mouse_df_filtered$conc_deposited_end_q25[i]
  vitro_internal_conc_q75    <- combined_mouse_df_filtered$conc_deposited_end_q75[i]
  vitro_internal_conc_q95    <- combined_mouse_df_filtered$conc_deposited_end_q95[i]
  
  density_liver = 1
  
  filtered_PK <- merged_PK[grepl(NP_core, merged_PK$NP_core), ]
  
  if (nrow(filtered_PK) > 0) {
    closest_row <- filtered_PK[which.min(abs(filtered_PK$Hydrodynamic_Size - HD)), ]
  } else {
    print(paste("No match found for NP_core:", NP_core, "instead find based on diameter"))
    closest_row <- merged_PK[which.min(abs(merged_PK$Hydrodynamic_Size - HD)), ]
  }
  
  # --- Extract DE quantiles ---
  DE_L_median <- closest_row$DE_Lt_id_g_median
  DE_L_q5     <- closest_row$DE_Lt_id_g_q5
  DE_L_q95    <- closest_row$DE_Lt_id_g_q95
  
  DE24_L_median <- closest_row$DE24_Lt_id_g_median
  DE24_L_q5     <- closest_row$DE24_Lt_id_g_q5
  DE24_L_q95    <- closest_row$DE24_Lt_id_g_q95
  
  liver_weight <- F_W_Liver * Body_weight * 1000
  n_mc <- 10000
  set.seed(42)
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
  # --- DE at max time ---
  result <- mc_ivive(
    bmd_median    = vitro_internal_conc_median,
    bmd_q5        = vitro_internal_conc_q5,
    bmd_q95       = vitro_internal_conc_q95,
    de_median     = DE_L_median,
    de_q5         = DE_L_q5,
    de_q95        = DE_L_q95,
    liver_weight  = liver_weight,
    density_liver = density_liver,
    body_weight   = Body_weight,
    n_mc          = n_mc
  )
  
  dose_q <- result$quantiles
  
  # --- DE at 24h ---
  result_24 <- mc_ivive(
    bmd_median    = vitro_internal_conc_median,
    bmd_q5        = vitro_internal_conc_q5,
    bmd_q95       = vitro_internal_conc_q95,
    de_median     = DE24_L_median,
    de_q5         = DE24_L_q5,
    de_q95        = DE24_L_q95,
    liver_weight  = liver_weight,
    density_liver = density_liver,
    body_weight   = Body_weight,
    n_mc          = n_mc
  )
  
  dose_24_q <- result_24$quantiles
  
  # --- Save results ---
  if (!any(is.na(dose_q))) {
    combined_mouse_df_filtered[i, "Estimated_Dose_q5_mg_kg_BMD"]     <- dose_q[1]
    combined_mouse_df_filtered[i, "Estimated_Dose_q25_mg_kg_BMD"]    <- dose_q[2]
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- dose_q[3]
    combined_mouse_df_filtered[i, "Estimated_Dose_q75_mg_kg_BMD"]    <- dose_q[4]
    combined_mouse_df_filtered[i, "Estimated_Dose_q95_mg_kg_BMD"]    <- dose_q[5]
    
    combined_mouse_df_filtered[i, "var_log_bmd10"]   <- result$var_bmd
    combined_mouse_df_filtered[i, "var_log_de"]      <- result$var_de
    combined_mouse_df_filtered[i, "pct_from_bmd10"]  <- result$pct_bmd
    combined_mouse_df_filtered[i, "pct_from_de"]     <- result$pct_de
  } else {
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- NA
    combined_mouse_df_filtered[i, "var_log_bmd10"]   <- result$var_bmd
    combined_mouse_df_filtered[i, "var_log_de"]      <- result$var_de
    combined_mouse_df_filtered[i, "pct_from_bmd10"]  <- result$pct_bmd
    combined_mouse_df_filtered[i, "pct_from_de"]     <- result$pct_de
  }
  
  if (!any(is.na(dose_24_q))) {
    combined_mouse_df_filtered[i, "Estimated_Dose_24_q5_mg_kg_BMD"]     <- dose_24_q[1]
    combined_mouse_df_filtered[i, "Estimated_Dose_24_q25_mg_kg_BMD"]    <- dose_24_q[2]
    combined_mouse_df_filtered[i, "Estimated_Dose_24_median_mg_kg_BMD"] <- dose_24_q[3]
    combined_mouse_df_filtered[i, "Estimated_Dose_24_q75_mg_kg_BMD"]    <- dose_24_q[4]
    combined_mouse_df_filtered[i, "Estimated_Dose_24_q95_mg_kg_BMD"]    <- dose_24_q[5]
    
    combined_mouse_df_filtered[i, "var_log_bmd10_24"]   <- result_24$var_bmd
    combined_mouse_df_filtered[i, "var_log_de_24"]      <- result_24$var_de
    combined_mouse_df_filtered[i, "pct_from_bmd10_24"]  <- result_24$pct_bmd
    combined_mouse_df_filtered[i, "pct_from_de_24"]     <- result_24$pct_de
  }
  
  combined_mouse_df_filtered[i, "closest_row_id"] <- closest_row$id
  #combined_mouse_df_filtered[i, "DE_L_id_g"]      <- DE_L_median
  #combined_mouse_df_filtered[i, "DE24_L_id_g"]    <- DE24_L_median
}

stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/"
save(combined_mouse_df_filtered, file = file.path(stanoutdir,"translated_dose_mouse_all.RData"))
write.csv(combined_mouse_df_filtered, file.path(stanoutdir,"translated_dose_mouse_all.csv"), row.names = FALSE)
