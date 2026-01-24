
library(dplyr)
library(ggplot2)
library(stringr)

folder = ".../mouse/bmd_results/"

#-------------0. setting-------
Body_weight = 0.02  # kg, Bodyweight (Davies 1993)
F_W_Liver  = 0.055  # (Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401)
F_W_Brain  = 0.017 # Unitless, Fractional brain tissue (Brown 1997)
F_W_Lung   = 0.007  # Unitless, Fractional lung tissue (Brown 1997)
F_W_Kidney = 0.017  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Spleen = 0.005  # Unitless®, Fractional spleen tissue (Brown 1997)


#----------1. PBPK model kinetic indicators ---------
predicted_paras = read.csv(".../pars_T_tot_mouse.csv")
preditced_PK = read.csv(".../PK_results_mouse.csv")
merged_PK <- merge(preditced_PK,predicted_paras , by.x = "id", by.y = "Folder", all = TRUE)


# ------- 2.1 load BMD results calculated from in vitro data directly -----

load("~.../combined_mouse_ec_filtered_results.RData")


# -------- 3. match the in vitro dosimetry model results with PROast results ---

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
    
  } else {
    print(paste("No match found for NP_core:", NP_core,"instead find based on diameter"))
    closest_row <- merged_PK[which.min(abs(merged_PK$Hydrodynamic_Size - HD)), ]
    # Extract the corresponding DE_L_id_g
    DE_L_id_g <- closest_row$DE_L_id_g
  }
  liver_weight = F_W_Liver * Body_weight * 1000 # gram as unit
  
  Dose_q5     <- (vitro_internal_conc_q5     * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q25    <- (vitro_internal_conc_q25    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_median <- (vitro_internal_conc_median * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q75    <- (vitro_internal_conc_q75    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q95    <- (vitro_internal_conc_q95    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  
  if (length(Dose_median) > 0) {
    print(1)
    combined_mouse_df_filtered[i, "Estimated_Dose_q5_mg_kg_BMD"]     <- Dose_q5     / (Body_weight * 1000)  # µg → mg/kg
    combined_mouse_df_filtered[i, "Estimated_Dose_q25_mg_kg_BMD"]    <- Dose_q25    / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- Dose_median / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_q75_mg_kg_BMD"]    <- Dose_q75    / (Body_weight * 1000)
    combined_mouse_df_filtered[i, "Estimated_Dose_q95_mg_kg_BMD"]    <- Dose_q95    / (Body_weight * 1000)
    } else {
    combined_mouse_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- NA
  }
  
}

stanoutdir = ".../mouse/"
save(combined_mouse_df_filtered, file = file.path(stanoutdir,"translated_dose_mouse_all.RData"))
