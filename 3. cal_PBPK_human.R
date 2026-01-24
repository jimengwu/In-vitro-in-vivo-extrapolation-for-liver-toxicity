library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' package
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
source("human_PBPK.R")
library(stringr)
library(tidyr)
library(tibble)
library(stringr)
library(readxl)

#-----------iv------
pred.mouse.iv <- function(pars,PDOSE,measurement_time,tstep) {
  
  ## Define the exposure scenario
  tinterval    = 1                                 ## hr, Time interval
  TDoses       = 1                                 ## Dose times, only one dose
  DOSE         = PDOSE*BW                          ## mg, amount of iv dose
  ex.iv<- ev(ID=1, amt= DOSE,                  ## Set up the exposure events
             ii=tinterval, addl=TDoses-1, 
             cmt="MBV", replicate = FALSE) 
  
  ## Set up the exposure time
  tsamp=tgrid(0,measurement_time,tstep)     ## Simulation time 24*7 hours (180 days)
  
  ## calculate the deposition volume
  out <- 
    mod %>% 
    param(pars) %>%
    update(atol=1e-50,maxsteps = 500000000) %>%
    mrgsim_d(data = ex.iv, tgrid=tsamp)
  
  ## save the calculated into data frame
  
  outdf = cbind.data.frame (Time       = out$time, 
                            CL=out$Liver_t,
                            CS = out$Spleen_t,
                            CK = out$Kidney_t,
                            Clung = out$Lung_t,
                            AUC_Kt_id_g     = out$AUC_Kt/(DOSE*10^6)*100,#from ng/g to %ID/g
                            AUC_Lt_id_g     = out$AUC_Lt/(DOSE*10^6)*100,
                            AUC_St_id_g     = out$AUC_St/(DOSE*10^6)*100,
                            AUC_Lut_id_g    = out$AUC_Lut/(DOSE*10^6)*100,
                            AUC_blood_id_g = out$AUC_blood/(DOSE*10^6)*100,
                            CLt_id_g        = out$Liver_t/(DOSE*10^6)*100,
                            CKt_id_g        = out$Kidney_t/(DOSE*10^6)*100,
                            CSt_id_g        = out$Spleen_t/(DOSE*10^6)*100,
                            CLungt_id_g     = out$Lung_t/(DOSE*10^6)*100,
                            CBlood_id_g      = out$Plasma/(DOSE*10^6)*100
  ) 
  return (list("outdf"  = outdf))
  
  return(out)
}


#-------------0. setting-------
BW = 73  # kg, Bodyweight (Davies 1993)

F_W_Liver  = 0.0257  # Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401
F_W_Brain  = 0.02    # Unitless, Fractional brain tissue (Brown 1997)
F_W_Lung   = 0.0076  # Unitless, Fractional lung tissue (Brown 1997)
F_W_Kidney = 0.0044  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Spleen = 0.0026  # Unitless®, Fractional spleen tissue (Brown 1997)

#----------- 1. calculate the human PBPK model ---------------

mod <- mcode ("human_PBPK", humanPBPK.code)

set.seed(5)

predicted_paras = read.csv(".../pars_T_tot.csv")
dataset_info <- read_excel(".../dataset_info.xlsx")

auc_results_list <- list()  # Initialize empty list

for (i in seq_len(nrow(predicted_paras))) {
  print(i)
  row <- predicted_paras[i, ]
  pars <- row[2:30]
  
  # Match case row from dataset_info
  case <- dataset_info %>% filter(id == row$Folder)
  
  # Extract dose and measurement time
  PDOSE <- case$Dose
  measurement_time <- case$Maximum_Measurement_Time
  
  # Set tstep conditionally
  tstep <- if (case$id == "FeO: Study2_41nm_4mg/kg") {
    0.5 / 60
  } else if (case$id %in% c(
    "GO: Study2_243nm_1mg/kg",
    "GO: Study2_914nm_1mg/kg_all",
    "GO: Study2_914nm_1mg/kg_w/o_CS"
  )) {
    1 / 60
  } else {
    min(1, measurement_time)
  }
  
  # Run the prediction
  R <- pred.mouse.iv(pars, PDOSE, measurement_time, tstep)
  
  # Extract AUC columns
  auc_cols <- grep("^AUC", names(R$outdf), value = TRUE)
  last_auc_row <- R$outdf[nrow(R$outdf), auc_cols, drop = FALSE]
  
  # Time-specific AUC extraction
  if (max(R$outdf$Time, na.rm = TRUE) >= 168) {
    auc24_row <- R$outdf %>%
      filter(Time == 24) %>%
      select(all_of(auc_cols)) %>%
      rename_with(~ str_replace_all(., "AUC", "AUC24"))
    
    auc168_row <- R$outdf %>%
      filter(Time == 168) %>%
      select(all_of(auc_cols)) %>%
      rename_with(~ str_replace_all(., "AUC", "AUC168"))
  } else {
    auc24_row <- data.frame(
      AUC24_Kt_id_g = NA, AUC24_Lt_id_g = NA, AUC24_St_id_g = NA,
      AUC24_Lut_id_g = NA, AUC24_blood_id_g = NA
    )
    auc168_row <- data.frame(
      AUC168_Kt_id_g = NA, AUC168_Lt_id_g = NA, AUC168_St_id_g = NA,
      AUC168_Lut_id_g = NA, AUC168_blood_id_g = NA
    )
  }
  
  # Max concentration values
  c_cols <- grep("^C", names(R$outdf), value = TRUE)
  max_c_values <- sapply(c_cols, function(col) {
    max_pos <- which.max(R$outdf[[col]])
    c(max = R$outdf[[col]][max_pos], time = R$outdf$Time[max_pos])
  })
  
  max_c_row <- as.data.frame(t(max_c_values)) %>%
    rownames_to_column(var = "id") %>%
    pivot_wider(names_from = id, values_from = -id)
  
  # Combine all extracted data
  result <- cbind(last_auc_row, auc24_row, auc168_row, max_c_row)
  result$time_last <- max(R$outdf$Time, na.rm = TRUE)

  # Combine with row info
  AUC_result <- cbind(row, result)
  
  # Save in list
  auc_results_list[[i]] <- AUC_result
}

# Combine all iterations into a single data frame
AUC_result_all <- bind_rows(auc_results_list)


# calculate the delivery efficiency based on the equation AUC/t_last
AUC_result_all$DE_L_id_g = AUC_result_all$AUC_Lt_id_g/AUC_result_all$time_last
AUC_result_all$DE_K_id_g = AUC_result_all$AUC_Kt_id_g/AUC_result_all$time_last
AUC_result_all$DE_S_id_g = AUC_result_all$AUC_St_id_g/AUC_result_all$time_last
AUC_result_all$DE_Lu_id_g = AUC_result_all$AUC_Lut_id_g/AUC_result_all$time_last

# caclualte the delivery efficiency based on the equation AUC/t_last
AUC_result_all$DE24_L_id_g = AUC_result_all$AUC24_Lt_id_g/24
AUC_result_all$DE24_K_id_g = AUC_result_all$AUC24_Kt_id_g/24
AUC_result_all$DE24_S_id_g = AUC_result_all$AUC24_St_id_g/24
AUC_result_all$DE24_Lu_id_g = AUC_result_all$AUC24_Lut_id_g/24

# caclualte the delivery efficiency based on the equation AUC/t_last
AUC_result_all$DE168_L_id_g = AUC_result_all$AUC168_Lt_id_g/168
AUC_result_all$DE168_K_id_g = AUC_result_all$AUC168_Kt_id_g/168
AUC_result_all$DE168_S_id_g = AUC_result_all$AUC168_St_id_g/168
AUC_result_all$DE168_Lu_id_g = AUC_result_all$AUC168_Lut_id_g/168


merged_AUC <- merge(dataset_info,AUC_result_all , by.x = "id", by.y = "Folder", all = TRUE)

write.csv(merged_AUC, ".../PK_results_human.csv", row.names = FALSE)

