
convert_HD_column <- function(df) {
  df %>%
    mutate(
      HD = na_if(HD, "-"),# Convert "-" to NA
      HD = as.numeric(HD)  # Convert to numeric (double)
    )
}

Exp_func <- function(x, a, b, c, d) {
  a * (c^(1 - exp(-(x / b)^d)))
}

Hill_func <- function (x,a,b,c,d){
  a*c^(x^d/(b^d+x^d))
}

get_dr_func <- function(modelname, a, b, c, d) {
  # Initialize dr_func and function_string as NULL (default values)
  dr_func <- NULL
  function_string <- NULL
  
  # Check for the model name and assign corresponding function and string
  if (grepl("H3", modelname) || grepl("H5", modelname)) {
    dr_func <- Hill_func  # Assign Hill function
    function_string <- paste0(
      "y = ", round(a,1), " * ", round(c,1), "^(x^", round(d,1), " / (", round(b,1), "^", round(d,1), " + x^", round(d,1), "))"
    )
    
  } else if (grepl("E3", modelname) || grepl("E5", modelname)) {
    dr_func <- Exp_func  # Assign Exp function
    function_string <- paste(
      "y = ",round(a,1), " * (", round(c,1), "^(1 - exp(-(x / ", round(b,1), ")^", round(d,1), ")))"
    )
    
    
  } else {
    stop("Model name not recognized!")
  }
  # Return the function and the string
  return(list(dr_func = dr_func, function_string = function_string))
  
}


plot_DR_curve <- function(x, y, a, b, c, d,BMD, modelname, 
                          datapoints_x, datapoints_y, xaxis_value="Concentration",
                          yaxis_value = "Response") {
  
  # Create the single_data dataframe for the scatter points
  single_data <- data.frame(
    x_value = datapoints_x,
    y_value = datapoints_y
  )
  
  # Create the data dataframe for the model curve
  line_data <- data.frame(
    x = x,
    y = y
  )
  
  
  # Call the function
  result <- get_dr_func(modelname, a, b, c, d)
  
  # Extract dr_func and function_string from the result
  dr_func <- result$dr_func
  function_string <- result$function_string
  
  
  # Create the plot
  p<- ggplot(line_data, aes(x = x, y = y)) +
    # Plot the model curve with a thicker line and color mapping for the legend
    geom_line(aes(color = "Model Curve"), linewidth = 1.5) +
    
    # Add scatter points from single_data with a different color for the legend
    geom_point(data = single_data, aes(x = x_value, y = y_value, color = "Data Points"), 
               size = 3, shape = 21, stroke = 0.5, fill = "#F1A09D") +
    
    # Add the vertical reference line for CED with a different color for the legend
    geom_segment(aes(x = BMD, y = 0, xend = BMD, yend = dr_func(BMD, a, b, c, d)), 
                 color = "black", linetype = "dashed", linewidth = 0.6) +
    geom_segment(aes(x = 0, y = dr_func(BMD, a, b, c, d), xend = BMD, yend = dr_func(BMD, a, b, c, d)), 
                 color = "black", linetype = "dashed", linewidth = 0.6) +
    
    # Add text annotation for CED value
    annotate("text", x = mean(line_data$x), y = 0.8*dr_func(BMD, a, b, c, d), 
             label = paste("CED =", round(BMD, 1)), 
             color = "black", size = 5, fontface = "bold") +
    
    # Add the equation as an annotation #todo 
    #annotate("text", x = BMD, y = 1.5 * dr_func(BMD, a, b, c, d), 
    #        label = function_string, 
    #        size = 5, color = "black") +
    
    
    # Logarithmic scales with nicely formatted axis ticks
    #scale_x_log10(expand = expansion(mult = c(0.05, 0.05))) + 
    #scale_y_log10() + 
    
    
    # Axis labels with bold text
    labs(x = xaxis_value, y = yaxis_value) +
    
    # Theme for a clean, publication-quality plot with grid lines
    theme_pubr(base_size = 16, base_family = "serif") + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),  # Add thick black border
      panel.grid.major = element_line(color = "gray70", linetype = "dashed", linewidth = 0.5),  # Major grid lines
      panel.grid.minor = element_line(color = "gray85", linetype = "dotted", linewidth = 0.4),  # Minor grid lines
      axis.text = element_text(size = 14, color = "black"),  # Large axis text
      axis.title = element_text(size = 16, face = "bold"),  # Bold axis labels
      legend.position = "top",  # Position legend at the top
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 14, color = "black"),  # Format legend text
      plot.background = element_blank()  # Remove any background shading
    ) +
    
    # Set custom colors for the legend
    scale_color_manual(values = c("Model Curve" = "#A1CDE2", "Data Points" = "black"))  # Set custom colors for legend
  return (p)
}




plot_vitro_toxicity <- function(data, y_col='Estimated_Dose_mg_kg',output_path = "plots/mouse/estimated_dose_vitro_mouse.pdf",plot_title="undeifined") {
  library(ggplot2)
  
  vitro_toxicity_plot <- ggplot(data, aes(x = NP, y = .data[[y_col]])) +
    geom_boxplot(fill = "#69b3a2", color = "#1b4f72", outlier.color = "#e74c3c", outlier.shape = 19) +
    scale_y_log10() +  # Apply log scale to the y-axis
    theme_minimal(base_size = 14) +  # Minimal theme with larger base font size
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"), # Rotate and bold x-axis labels
      axis.text.y = element_text(size = 12, face = "bold"),                        # Bold y-axis labels
      axis.title.x = element_text(size = 14, face = "bold"),                       # Bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),                       # Bold y-axis title
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),            # Centered bold plot title
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),      # Light dashed grid lines
      panel.grid.minor = element_blank(),                                           # Remove minor grid lines
      panel.background = element_rect(fill = "white", color = "black", size = 1.5) # Add black border around plot area
    ) +
    labs(
      title = plot_title,
      x = "Substance",
      y = "Estimated Dose (mg/kg)"
    )
  
  # Save the plot
  ggsave(output_path, plot = vitro_toxicity_plot, width = 10, height = 7)
  
  return(vitro_toxicity_plot) # Return the plot object
}



plot_vitro_toxicity_category <- function(data, 
                                         y_col='Estimated_Dose_mg_kg',
                                         output_path = "plots/mouse/estimated_dose_vitro_mouse.pdf",
                                         plot_title) {
  library(ggplot2)
  
  vitro_toxicity_plot <- ggplot(data, aes(x = NP, y = .data[[y_col]])) +
    geom_boxplot(aes(fill = assay_category),  # Directly use assay_category for legend mapping
                 color = "#1b4f72", outlier.color = "#e74c3c", outlier.shape = 19) +
    scale_y_log10() +  # Log scale for better visualization
    scale_fill_manual(
      values = c("viability" = "#69b3a2", "oxidative_stress" = "#ffa07a"),  # Define two colors
      name = "Assay Category"  # Legend title
    ) +
    theme_minimal(base_size = 14) +  # Minimal theme with a readable font size
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  # Rotate and bold x-axis labels
      axis.text.y = element_text(size = 12, face = "bold"),  # Bold y-axis labels
      axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),  # Bold y-axis title
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered bold plot title
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),  # Light dashed grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white", color = "black", size = 1.5)  # Add black border around plot area
    ) +
    labs(
      title = plot_title,
      x = "Substance",
      y = "Estimated Dose (mg/kg)"
    )
  
  # Save the plot
  ggsave(output_path, plot = vitro_toxicity_plot, width = 10, height = 7)
  
  return(vitro_toxicity_plot) # Return the plot object
}


plot_mouse_comparison <- function(data, output_path = "plots/mouse/mouse_comparison_plot.pdf",plot_title) {
  library(ggplot2)
  
  mouse_comparison_plot <- ggplot(data, aes(x = Substance, y = Dose, fill = Source)) +
    geom_boxplot(
      position = position_dodge(width = 0.8),  # Dodge for side-by-side positioning
      color = "black", outlier.color = "red", outlier.shape = 19
    ) +
    scale_y_log10() +  # Apply log scale to the y-axis
    scale_fill_manual(values = c("in vitro" = "#69b3a2", "in vivo" = "#ffa07a")) +  # Custom colors
    theme_minimal(base_size = 14) +  # Minimal theme with larger base font size
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  # Rotate and bold x-axis labels
      axis.text.y = element_text(size = 12, face = "bold"),                         # Bold y-axis labels
      axis.title.x = element_text(size = 14, face = "bold"),                        # Bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),                        # Bold y-axis title
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),             # Centered bold plot title
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),       # Light dashed grid lines
      panel.grid.minor = element_blank(),                                           # Remove minor grid lines
      legend.title = element_blank(),                                                # Remove legend title
      panel.background = element_rect(fill = "white", color = "black", size = 1.5) # Add black border around plot area
    ) +
    labs(
      title=plot_title,
      x = "Substance",
      y = "Dose (mg/kg)"
    )
  
  # Save the plot
  ggsave(output_path, plot = mouse_comparison_plot, width = 10, height = 7)
  
  return(mouse_comparison_plot) # Return the plot object
}


Pred_auc <- function (pars,PDOSE,tstep){
  
  
  BW           = 73                              ## kg, body weight
  tinterval    = 1                                 ## hr, Time interval for input
  TDoses       = 1                                 ## Dose times, only one dose
  
  ## Get out of log domain
  # pars is the entire parameter set; pars [-which_sig] means to keep 
  # parameters without "sig" only, and then do exp transformation, then reassign to pars
  pars <- lapply(pars,exp) 
  
  ## Repeat dose exposure scenario: 
  DOSE    = PDOSE*BW            ## mg; amount of oral dose
  
  ex         <- ev(ID=1, amt= DOSE, ii=tinterval, 
                   addl=TDoses-1, cmt="MBV", replicate = FALSE)
  
  
  ## set up the exposure time
  ## Simulated for 24*365 hours after dosing, but only obtained data at 24 h
  tsamp     = tgrid(0,max(Obs.df$Time),tstep)          
  
  
  ## Get a prediction
  # The code can produce time-dependent NSC values, but at time = 0, 
  # NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  out <- 
    mod %>%
    param(pars) %>%
    update(atol = 1E-80,maxsteps = 5000000)%>%
    mrgsim_d(data = ex, tgrid = tsamp)%>%
    filter(time!=0) 
  
  outdf = cbind.data.frame (Time       = out$time, 
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
  
}


generate_combined_results <- function(list_object, results_df) {
  ex_info <- data.frame(
    Paper_ID = sapply(list_object, function(x) unique(x$`Study_ID`)),
    NP = sapply(list_object, function(x) unique(x$`Substance name`)),
    NP_core = sapply(list_object, function(x) paste(unique(x$`Substance core`), collapse = ", ")),
    #HD = sapply(list_object, function(x) unique(x$`Hydrodynamic size (nm)`)),
    assay =sapply(list_object, function(x) unique(x$`Tested assay`)), 
    HD = sapply(list_object, function(x) paste(unique(x$`Hydrodynamic size (nm)`), collapse = ", ")),
    exposure_time = sapply(list_object, function(x) unique(x$`In vitro exposure time[h]`)),
    Shape = sapply(list_object, function(x) unique(x$`Shape`)),
    Average_size_nm = sapply(list_object, function(x) unique(x$`Average size (nm)`)),
    Z_potential_mV = sapply(list_object, function(x) unique(x$`Z potential [mV]`)[1]),
    #Z_potential_mV = sapply(list_object, function(x) paste(unique(x$`Z potential [mV]`), collapse = ", ")),
    Cell_type = sapply(list_object, function(x) unique(x$`Cell type`)),
    media_type = sapply(list_object, function(x) unique(x$`In vitro media`)),
    stringsAsFactors = FALSE
  )
  # Add sequential study_id
  ex_info$study <- seq_len(nrow(ex_info))
  
  #combined_df <- cbind(ex_info, results_df)
  combined_df <-merge(ex_info, results_df, by = "study")
  return(combined_df)
}



plot_ecx_human <- function(df,
                           ec_prefix = "ec10",   
                           np_col    = "NP",
                           cell_col  = "Cell_type",
                           exp_col   = "general_experiment_type") {
  
  # build column names for the EC summaries
  ec_med  <- paste0(ec_prefix, "_median")
  ec_q5   <- paste0(ec_prefix, "_q5")
  ec_q25  <- paste0(ec_prefix, "_q25")
  ec_q75  <- paste0(ec_prefix, "_q75")
  ec_q95  <- paste0(ec_prefix, "_q95")
  
  
  # Prepare plotting data with ordered combo_label
  df_plot <- df %>%
    mutate(
      combo_label = paste(.data[[np_col]], .data[[exp_col]], .data[[cell_col]], sep = " | "),
      combo_label = fct_reorder(combo_label, .data[[ec_med]], .fun = min, .na_rm = TRUE)
    )
  
  df_plot <- df %>%
    mutate(
      combo_label = paste(.data[[np_col]], .data[[cell_col]], sep = " | "),
      combo_label = fct_reorder(combo_label, .data[[ec_med]], .fun = min, .na_rm = TRUE)
    )
  
  ggplot(df_plot, aes(y = combo_label)) +
    # 90% interval (q5–q95) – context
    geom_errorbarh(aes(xmin = .data[[ec_q5]], xmax = .data[[ec_q95]]),
                   height = 0.2, color = "#808080", size = 1,alpha=1) +
    # 50% interval (q25–q75) – emphasis
    #geom_errorbarh(aes(xmin = .data[[ec_q25]], xmax = .data[[ec_q75]]),
    #               height = 0.2, color = "#96cccb", size = 3) + #4198AC
    geom_segment(
      aes(x = .data[[ec_q25]], 
          xend = .data[[ec_q75]], 
          y = combo_label, 
          yend = combo_label),
      #color = "#96cccb", fill="#e66d50"
      color="#808080",
      alpha = 0.7,
      size = 3
    )+
    #geom_point(aes(x = .data[[ec_med]], shape = .data[[exp_col]]),
    #           color = "black", fill="#BBAADE", size = 2, stroke = 1.5,alpha = 0.95) +
    # scale_shape_manual(values = c(21, 22, 23, 24),
    #                  name   = "Experiment type") +
    # Median point
    geom_point(aes(x = .data[[ec_med]],  fill = .data[[exp_col]]), color="black",
               shape = 21,  size = 3, alpha = 0.95,stroke=1.5) +
    #scale_fill_manual(values = c("#3E72B6", "#BF5A2D", "#6C3687", "#DAB039")) +
    scale_fill_manual(name = "Experiment type",values = c( "#BF5A2D","#8BAFCD", "#6C3687","#DAB039")) +
    scale_x_log10(labels = label_scientific(digits = 1.8),limit = c(1e-5,1e4)) +
    labs(
      #x = bquote(bold(EC[.(gsub("ec", "", ec_prefix))]) ~ "(" * mu * "g/ml)"),
      x = bquote(bolditalic("In vitro") ~ bold(BMD[.(gsub("ec", "", ec_prefix))] ~ "(" * mu * "g/ml)")),
      y = "NP core | cell type",
      color = "Experiment type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y  = element_text(size = 11, color = "black"),
      axis.text.x  = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_line(color = "gray85", linetype = "dashed"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor   = element_blank(),
      
      legend.position      = c(0.01, 0.99),      # bottom-left
      legend.justification = c(0, 1),            # anchor legend box to bottom-left
      legend.title         = element_text(size = 11, face = "bold"),
      legend.text          = element_text(size = 10),
      legend.box.background = element_rect(color = "black", size = 0.7, linetype = "solid"),
      
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background   = element_rect(fill = "white", color = NA)
    )
  
  
  
}



shorten_assay_names <- function(df, assay_col = "assay") {
  df %>%
    mutate(
      assay_short = case_when(
        str_detect(.data[[assay_col]], regex("alanine aminotransferase.*ALT", ignore_case = TRUE)) ~ "ALT",
        str_detect(.data[[assay_col]], regex("alkaline phosphatase.*ALP", ignore_case = TRUE)) ~ "ALP",
        str_detect(.data[[assay_col]], regex("aspartate transaminase.*AST", ignore_case = TRUE)) ~ "AST",
        str_detect(.data[[assay_col]], regex("microalblmin.*MIA", ignore_case = TRUE)) ~ "MIA",
        str_detect(.data[[assay_col]], regex("protein assay", ignore_case = TRUE)) ~ "Protein",
        str_detect(.data[[assay_col]], regex("total bilirubin.*T BIL", ignore_case = TRUE)) ~ "T BIL",
        str_detect(.data[[assay_col]], regex("triglycerides.*TRIG", ignore_case = TRUE)) ~ "TRIG",
        str_detect(.data[[assay_col]], regex("catalase.*CAT", ignore_case = TRUE)) ~ "CAT",
        str_detect(.data[[assay_col]], regex("gamma glutamyltranspeptidase.*GGT", ignore_case = TRUE)) ~ "GGT",
        str_detect(.data[[assay_col]], regex("glucose 6-phosphate dehydrogenase.*G6PDH", ignore_case = TRUE)) ~ "G6PDH",
        str_detect(.data[[assay_col]], regex("glutathione peroxidase.*GPx", ignore_case = TRUE)) ~ "GPx",
        str_detect(.data[[assay_col]], regex("glutathione reductase.*GRD", ignore_case = TRUE)) ~ "GRD",
        str_detect(.data[[assay_col]], regex("reduced glutathione.*GSH", ignore_case = TRUE)) ~ "GSH",
        str_detect(.data[[assay_col]], regex("thioredoxin reductase.*thrr", ignore_case = TRUE)) ~ "THRR",
        str_detect(.data[[assay_col]], regex("malondialdehyde.*MDA", ignore_case = TRUE)) ~ "MDA",
        str_detect(.data[[assay_col]], regex("protein carbonyl.*PC", ignore_case = TRUE)) ~ "PC",
        str_detect(.data[[assay_col]], regex("ROS.*DCFH-DA|% of ROS|fluorescence|DCFDA", ignore_case = TRUE)) ~ "ROS",
        str_detect(.data[[assay_col]], regex("superoxide dismutase.*SOD", ignore_case = TRUE)) ~ "SOD",
        TRUE ~ .data[[assay_col]]  # fallback: retain full string
      )
    )
}

shorten_cell_types <- function(df, cell_col = "Cell_type") {
  df %>%
    mutate(
      !!cell_col := case_when(
        str_detect(.data[[cell_col]], regex("human hepatoblastoma C3A cell line", ignore_case = TRUE)) ~ "C3A",
        str_detect(.data[[cell_col]], regex("primary human hepatocytes", ignore_case = TRUE)) ~ "PHH",
        str_detect(.data[[cell_col]], regex("Hepa1-6", ignore_case = TRUE)) ~ "Hepa 1-6",
        str_detect(.data[[cell_col]], regex("Hepa 1-6", ignore_case = TRUE)) ~ "Hepa 1-6",
        TRUE ~ .data[[cell_col]]  # fallback: retain full string
      )
    )
}

library(tidyverse)
library(scales)

# =============================================================================
# MAIN TEXT — Option 5: Strip + shaded IQR band
# Shows individual experiment points jittered over a shaded IQR band per group
# =============================================================================

plot_ecx_human_main <- function(df,
                                ec_prefix = "ec10",
                                np_col    = "NP",
                                cell_col  = "Cell_type",
                                exp_col   = "general_experiment_type",
                                min_n     = 3) {
  
  # Column names
  ec_med  <- paste0(ec_prefix, "_median")
  ec_q5   <- paste0(ec_prefix, "_q5")
  ec_q25  <- paste0(ec_prefix, "_q25")
  ec_q75  <- paste0(ec_prefix, "_q75")
  ec_q95  <- paste0(ec_prefix, "_q95")
  
  # --- Build combo label (NP | cell type) and order by group median ---
  # Ordering uses ALL groups, regardless of n
  df_plot <- df %>%
    mutate(combo_label = paste(.data[[np_col]], .data[[cell_col]], sep = " | ")) %>%
    mutate(
      combo_label = fct_reorder(combo_label, .data[[ec_med]],
                                .fun = median, .na_rm = TRUE)
    )
  
  # --- Subset used for boxes only: groups with n >= min_n non-missing values ---
  df_box <- df_plot %>%
    group_by(combo_label) %>%
    filter(sum(!is.na(.data[[ec_med]])) >= min_n) %>%
    ungroup()
  
  ggplot() +
    # Layer 1: Box-whisker plot, only for groups with n >= min_n
    geom_boxplot(
      data = df_box,
      aes(x = .data[[ec_med]], y = combo_label),
      fill = "#F0EBE1", alpha = 1, color = "gray40",
      outlier.shape = NA, width = 0.5, linewidth = 0.4
    ) +
    # Layer 2: Individual experiment medians for ALL groups
    geom_point(
      data = df_plot,
      aes(x = .data[[ec_med]], y = combo_label,
          fill = .data[[exp_col]]),
      shape = 21, size = 2.5, stroke = 0.8, color = "black",
      position = position_jitter(width = 0, height = 0.15, seed = 42)
    ) +
    scale_fill_manual(name = "Experiment type",
                      values = c("Cytokine"              = "#BF5A2D",
                                 "Cytotoxicity"          = "#8BAFCD",
                                 "Genotoxicity"          = "#6C3687",
                                 "Metabolic Cell Stress" = "#DAB039")) +
    scale_y_discrete(drop = FALSE) +
    scale_x_log10(
      labels = label_scientific(digits = 1),
      limits = c(1e-3, 1e4),
      breaks = 10^(-3:4)
    ) +
    labs(
      x = bquote(bolditalic("In vitro") ~
                   bold(BMD[.(gsub("ec", "", ec_prefix))] ~ "(" * mu * "g/mL)")),
      y = "NP core | cell type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y    = element_text(size = 11, color = "black"),
      axis.text.x    = element_text(size = 11, color = "black"),
      axis.title.x   = element_text(size = 12, face = "bold"),
      axis.title.y   = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_line(color = "gray85", linetype = "dashed"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position      = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.title         = element_text(size = 11, face = "bold"),
      legend.text          = element_text(size = 10),
      legend.box.background = element_rect(color = "black", linewidth = 0.5),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white", color = NA)
    )
}


# =============================================================================
# SUPPORTING INFORMATION — Option 2: Faceted interval plot
# One panel per nanoparticle, each experiment on its own row with full CI
# =============================================================================

plot_ecx_human_si <- function(df,
                              ec_prefix = "ec10",
                              np_col    = "NP",
                              cell_col  = "Cell_type",
                              exp_col   = "general_experiment_type",
                              id_col    = NULL,
                              size_col  = NULL,
                              time_col  = NULL,
                              study_col = "study") {
  
  # Column names
  ec_med  <- paste0(ec_prefix, "_median")
  ec_q5   <- paste0(ec_prefix, "_q5")
  ec_q25  <- paste0(ec_prefix, "_q25")
  ec_q75  <- paste0(ec_prefix, "_q75")
  ec_q95  <- paste0(ec_prefix, "_q95")
  
  # --- Build facet label (NP | cell type) ---
  df_plot <- df %>%
    mutate(
      facet_label = paste(.data[[np_col]], .data[[cell_col]], sep = " | ")
    )

  # --- Build row label with remaining metadata (size, time, study) ---
  df_plot <- df_plot %>%
    mutate(row_label = paste0(.data[["experiment_type"]], " No.", .data[[study_col]]))
  
  # Append particle size if column exists
  if (!is.null(size_col) && size_col %in% names(df)) {
    df_plot <- df_plot %>%
      mutate(row_label = paste0(row_label, " | ", round(.data[[size_col]]), " nm"))
  }
  
  # Append exposure time if column exists
  if (!is.null(time_col) && time_col %in% names(df)) {
    df_plot <- df_plot %>%
      mutate(row_label = paste0(row_label, " | ", .data[[time_col]], "h"))
  }
  
  # Append experiment ID if available (for uniqueness)
  if (!is.null(id_col) && id_col %in% names(df)) {
    df_plot <- df_plot %>%
      mutate(row_label = paste0(row_label, " (", .data[[id_col]], ")"))
  }
  
  # Deduplicate labels if needed and order by median BMD10
  df_plot <- df_plot %>%
    mutate(
      row_label = make.unique(row_label, sep = " "),
      row_label = fct_reorder(row_label, .data[[ec_med]], .na_rm = TRUE)
    )
  
  # --- Compute summary stats per facet panel for background boxplot ---
  facet_summary <- df_plot %>%
    group_by(facet_label) %>%
    summarise(
      box_med   = median(.data[[ec_med]], na.rm = TRUE),
      box_q25   = quantile(.data[[ec_med]], 0.25, na.rm = TRUE),
      box_q75   = quantile(.data[[ec_med]], 0.75, na.rm = TRUE),
      box_lower = max(min(.data[[ec_med]], na.rm = TRUE),
                      quantile(.data[[ec_med]], 0.25, na.rm = TRUE) -
                        1.5 * IQR(.data[[ec_med]], na.rm = TRUE)),
      box_upper = min(max(.data[[ec_med]], na.rm = TRUE),
                      quantile(.data[[ec_med]], 0.75, na.rm = TRUE) +
                        1.5 * IQR(.data[[ec_med]], na.rm = TRUE)),
      .groups = "drop"
    )
  
  ggplot(df_plot, aes(y = row_label)) +
    geom_vline(
      data = facet_summary,
      aes(xintercept = box_lower),
      color = "gray70", linewidth = 0.3, linetype = "dotted"
    ) +
    geom_vline(
      data = facet_summary,
      aes(xintercept = box_upper),
      color = "gray70", linewidth = 0.3, linetype = "dotted"
    ) +
    geom_boxplot(
      aes(y      = row_label,
          xmin    = .data[[ec_q5]],
          xlower  = .data[[ec_q25]],
          xmiddle = .data[[ec_med]],
          xupper  = .data[[ec_q75]],
          xmax    = .data[[ec_q95]],
          fill    = .data[[exp_col]]),
      stat = "identity",
      color = "gray40",
      width = 0.6, linewidth = 0.4
    ) +
    facet_wrap(
      ~ facet_label,
      scales = "free_y",
      ncol   = 2
    ) +
    scale_fill_manual(name = "Experiment type",
                      values = c( "Cytokine"= "#BF5A2D",
                                            "Cytotoxicity"= "#8BAFCD", 
                                            "Genotoxicity"="#6C3687",
                                            "Metabolic Cell Stress"="#DAB039")) +
    scale_x_log10(
      labels = label_scientific(digits = 1),
      limits = c(1e-3, 1e4),
      breaks = 10^(-3:4)
    ) +
    labs(
      x = bquote(bolditalic("In vitro") ~
                   bold(BMD[.(gsub("ec", "", ec_prefix))] ~ "(" * mu * "g/mL)")),
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y    = element_text(size = 8, color = "black"),
      axis.text.x    = element_text(size = 9, color = "black"),
      axis.title.x   = element_text(size = 11, face = "bold"),
      strip.text      = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray70",
                                      linewidth = 0.5),
      panel.grid.major.x = element_line(color = "gray85", linetype = "dashed"),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      legend.position = "bottom",
      legend.title    = element_text(size = 10, face = "bold"),
      legend.text     = element_text(size = 9),
      panel.border     = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing    = unit(1, "lines")
    )
}

