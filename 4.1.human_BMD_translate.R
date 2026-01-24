# this function runs the PBPK model get the deposition fraction after 24 hours 
# and then combine it together with the in vitro deposition fraction got from python 
# try to get the deposited amount before calculating the benchmark dose, instead, we calculate the 
# human equivalent dose first and then get the bmd calculation
library(stringr)
library(dplyr)
library(ggplot2)
source("with_stan/helper_function.R")

folder = ".../bmd_results/"

#-------------0. setting-------
Body_weight = 73  # kg, Bodyweight (Davies 1993)
F_W_Liver  = 0.0257  # (Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401)
F_W_Brain  = 0.02 # Unitless, Fractional brain tissue (Brown 1997)
F_W_Lung   = 0.0076  # Unitless, Fractional lung tissue (Brown 1997)
F_W_Kidney = 0.0044  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Spleen = 0.0026  # Unitless®, Fractional spleen tissue (Brown 1997)


#----------1. PBPK model kinetic indicators ---------
preditced_human_PK = read.csv(".../PK_results.csv")
# load the processed effective concentration value
load("results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData")

combined_human_df_filtered$Fraction_deposited_end = 1
combined_human_df_filtered$conc_deposited_end_median <- 
  combined_human_df_filtered$ec10_median * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q5  <- combined_human_df_filtered$ec10_q5  * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q25 <- combined_human_df_filtered$ec10_q25 * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q75 <- combined_human_df_filtered$ec10_q75 * combined_human_df_filtered$Fraction_deposited_end
combined_human_df_filtered$conc_deposited_end_q95 <- combined_human_df_filtered$ec10_q95 * combined_human_df_filtered$Fraction_deposited_end


#--------------2. merge the in vitro model results_vitro with the PBPK model -------

#The sinusoidal cells are the predominant non-parenchymal cells, comprising about 35% of the total cell number and about 17% of the total volume of the liver.
#These cells are divided into sinusoidal endothelial cells (44%), Kupffer cells (33%), stellate cells (10-25%) and hepatic NK cells (5%)

library(dplyr)

combined_human_df_filtered <- combined_human_df_filtered %>%
  mutate(HD = ifelse(is.na(HD), Average_size_nm, HD))

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
  
  
  # Filter AUC_result_all to only rows with matching NP_core #todo: ZnO is not found,
  filtered_PK <- preditced_human_PK[grepl(NP_core, preditced_human_PK$NP_core), ]
  
  if (nrow(filtered_PK) > 0) {
    # Find the row with the closest Hydrodynamic_Size
    closest_row <- preditced_human_PK[which.min(abs(preditced_human_PK$Hydrodynamic_Size - HD)), ]
    # Extract the corresponding DE_L_id_g
    DE_L_id_g <- closest_row$DE_L_id_g
    DE24_L_id_g <- closest_row$DE24_L_id_g
    
  } else {
    print(paste("No match found for NP_core:", NP_core,"instead find based on diameter"))
    closest_row <- preditced_human_PK[which.min(abs(preditced_human_PK$Hydrodynamic_Size - HD)), ]
    # Extract the corresponding DE_L_id_g
    DE_L_id_g <- closest_row$DE_L_id_g
    DE24_L_id_g <- closest_row$DE24_L_id_g
  }
  
  liver_weight = F_W_Liver * Body_weight * 1000 # gram as unit
  # ug/cm3 for vitro_internal_conc
  
  Dose_q5     <- (vitro_internal_conc_q5     * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q25    <- (vitro_internal_conc_q25    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_median <- (vitro_internal_conc_median * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q75    <- (vitro_internal_conc_q75    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  Dose_q95    <- (vitro_internal_conc_q95    * liver_weight / density_liver) / (DE_L_id_g * liver_weight / 100)
  
  Dose_24_q5     <- (vitro_internal_conc_q5     * liver_weight / density_liver) / (DE24_L_id_g * liver_weight / 100)
  Dose_24_q25    <- (vitro_internal_conc_q25    * liver_weight / density_liver) / (DE24_L_id_g * liver_weight / 100)
  Dose_24_median <- (vitro_internal_conc_median * liver_weight / density_liver) / (DE24_L_id_g * liver_weight / 100)
  Dose_24_q75    <- (vitro_internal_conc_q75    * liver_weight / density_liver) / (DE24_L_id_g * liver_weight / 100)
  Dose_24_q95    <- (vitro_internal_conc_q95    * liver_weight / density_liver) / (DE24_L_id_g * liver_weight / 100)
  
  if (length(Dose_median) > 0) {
    combined_human_df_filtered[i, "Estimated_Dose_q5_mg_kg_BMD"]     <- Dose_q5     / (Body_weight * 1000)  # µg → mg/kg
    combined_human_df_filtered[i, "Estimated_Dose_q25_mg_kg_BMD"]    <- Dose_q25    / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_median_mg_kg_BMD"] <- Dose_median / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_q75_mg_kg_BMD"]    <- Dose_q75    / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_q95_mg_kg_BMD"]    <- Dose_q95    / (Body_weight * 1000)
    
    combined_human_df_filtered[i, "Estimated_Dose_24_q5_mg_kg_BMD"]     <- Dose_24_q5     / (Body_weight * 1000)  # µg → mg/kg
    combined_human_df_filtered[i, "Estimated_Dose_24_q25_mg_kg_BMD"]    <- Dose_24_q25    / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_24_median_mg_kg_BMD"] <- Dose_24_median / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_24_q75_mg_kg_BMD"]    <- Dose_24_q75    / (Body_weight * 1000)
    combined_human_df_filtered[i, "Estimated_Dose_24_q95_mg_kg_BMD"]    <- Dose_24_q95    / (Body_weight * 1000)
    
    
  } else {
    combined_human_df_filtered[i, "Estimated_Dose_mg_kg_BMD"] <- NA
  }
  
}

