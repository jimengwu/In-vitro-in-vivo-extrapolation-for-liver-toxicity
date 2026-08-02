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
pred.iv.mouse <- function(mod,pars,PDOSE,measurement_time,tstep) {
  
  ## Get out of log domain
  F_W_Liver  = 0.055  
  F_W_Brain  = 0.017 
  F_W_Lung   = 0.007  
  F_W_Kidney = 0.017
  F_W_Spleen = 0.005 
  F_W_Plasma = 0.029  
  
  BW           = 0.02                              ## kg, body weight
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


pred.iv.human <- function(mod,pars,PDOSE,measurement_time,tstep) {
  
  F_W_Liver  = 0.0257  # Unitless, Fractional liver tissue (Brown 1997) table 4, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401
  F_W_Brain  = 0.02    # Unitless, Fractional brain tissue (Brown 1997)
  F_W_Lung   = 0.0076  # Unitless, Fractional lung tissue (Brown 1997)
  F_W_Kidney = 0.0044  # Unitless, Fractional kidney tissue (Brown 1997)
  F_W_Spleen = 0.0026  # Unitless®, Fractional spleen tissue (Brown 1997)
  
  BW           = 73                              ## kg, body weight
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

#----------- 1. calculate the human PBPK model ---------------

mod_human <- mcode ("human_PBPK", humanPBPK.code)
mod_mouse <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)

#R <- pred.mouse.iv(pars, PDOSE, measurement_time, tstep)

ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg","GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")

predicted_paras = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/pars_T_tot.csv")
dataset_info <- read_excel("~/work/code/Mouse-general-PBPK/dataset/tk/mouse/dataset_info.xlsx")
mouse_AUC = read.csv("~/work/code/Mouse-general-PBPK/plots/auc_results.csv")
human_results <- data.frame()
for (i in 1:length(ls_np_name)) {
  row  <- predicted_paras[i, ]
  case <- dataset_info %>% filter(id == row$Folder)
  PDOSE <- case$Dose
  measurement_time <- case$Maximum_Measurement_Time
  
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
  
  # --- Run both models with bestpar ---
  pars_mouse_best <- as.numeric(row[2:30])
  names(pars_mouse_best) <- colnames(row)[2:30]
  names(pars_mouse_best) <- gsub("K_max","K_uptake", names(pars_mouse_best))
  
  R_mouse <- tryCatch(pred.iv.mouse(mod_mouse,pars_mouse_best, PDOSE, measurement_time, tstep),
                      error = function(e) NULL)
  R_human <- tryCatch(pred.iv.human(mod_human,pars_mouse_best, PDOSE, measurement_time, tstep),
                      error = function(e) NULL)
  mouse_de24 <- R_mouse$outdf[which.min(abs(R_mouse$outdf$Time - 24)), "AUC_Lt_id_g"] / 24
  human_de24 <- R_human$outdf[which.min(abs(R_human$outdf$Time - 24)), "AUC_Lt_id_g"] / 24
  ratio_24 <- human_de24 / mouse_de24
  
  mouse_de <- R_mouse$outdf[which.min(abs(R_mouse$outdf$Time - measurement_time)), "AUC_Lt_id_g"] / 24
  human_de <- R_human$outdf[which.min(abs(R_human$outdf$Time - measurement_time)), "AUC_Lt_id_g"] / 24
  ratio_max <- human_de / mouse_de
  
  mouse_row <- mouse_AUC %>% filter(id == row$Folder)
  if (nrow(mouse_row) == 0) {
    cat("  WARNING: No mouse MCMC result for", row$Folder, "- skipping\n")
    next
  }

  result_row <- data.frame(
    id                  = row$Folder,
    ratio_liver_24      = ratio_24,
    human_DE24_L_q5    = mouse_row$DE24_Lt_id_g_q5    * ratio_24,
    human_DE24_L_q25    = mouse_row$DE24_Lt_id_g_q25    * ratio_24,
    human_DE24_L_median = mouse_row$DE24_Lt_id_g_median * ratio_24,
    human_DE24_L_q75    = mouse_row$DE24_Lt_id_g_q75    * ratio_24,
    human_DE24_L_q95    = mouse_row$DE24_Lt_id_g_q95    * ratio_24,
    
    human_DE_L_q5    = mouse_row$DE_Lt_id_g_q5    * ratio_max,
    human_DE_L_q25    = mouse_row$DE_Lt_id_g_q25    * ratio_max,
    human_DE_L_median = mouse_row$DE_Lt_id_g_median * ratio_max,
    human_DE_L_q75    = mouse_row$DE_Lt_id_g_q75    * ratio_max,
    human_DE_L_q95    = mouse_row$DE_Lt_id_g_q95    * ratio_max
  )
  
  human_results <- rbind(human_results, result_row)
  cat("  Ratio:", round(ratio_24, 4), "\n")
}

human_results
human_results_merged <- human_results %>%
  left_join(dataset_info, by = c("id" = "id"))

write.csv(human_results_merged, "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/PK_results_ratio.csv", row.names = FALSE)

