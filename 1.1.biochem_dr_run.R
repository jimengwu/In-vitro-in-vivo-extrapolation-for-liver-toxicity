library(cmdstanr)
library(ggplot2)
library(dplyr)
library(purrr)
library(posterior)

capitalize_if_needed <- function(x) {
  ifelse(grepl("^[a-z]", x), 
         paste0(toupper(substring(x, 1, 1)), substring(x, 2)), 
         x)
}


plot_fit_initial <- function(x_obs,y_obs){
  hill_function_bc <- function(x, a, b,c, d) {
    a *(1 +(c-1)*(x^d / (b^d + x^d)))
  }
  
  neg_log_lik <- function(params, x, y) {
    a <- params[1]; b <- params[2]; c <- params[3];d <- params[4]; sigma <- params[5]
    y_hat <- hill_function_bc(x, a, b, c, d)
    -sum(dnorm(y, mean = y_hat, sd = sigma, log = TRUE))
  }
  
  runs <- map(1:5000, ~ {
    init <- c(
      a = runif(1, min(y_obs)/100, max(y_obs)*100), # baseline response
      # b = runif(1, 0.001, 100), # hill-ced model 
      b = runif(1, min(x_obs[x_obs > 0])/10, max(x_obs)*10), # for H5 model only,  inflection point (EC50 or BMD-like — dose at which half-max effect occurs)
      c = runif(1, 1, (max(y_obs)/min(y_obs))*5), # response scaling factor — asymptote at high dose, relative to baseline
      d = runif(1, 0.01, 10),  # slope or Hill coefficient — steepness of curve
      sigma = runif(1, 0.01, 0.1))
    optim(par = init, fn = neg_log_lik, x = x_obs, y = y_obs, method = "L-BFGS-B", 
          lower = c(min(y_obs)/100, min(x_obs[x_obs > 0])/10, 1, 0.01,1e-6), 
          upper = c(max(y_obs)*100, max(x_obs)*10, (max(y_obs)/min(y_obs))*5, 10, 1))
  })
  
  
  fit_opt <- runs[[which.min(purrr::map_dbl(runs, ~ .x$value))]]
  
  plot = plot_obs_pred(x_obs,y_obs,fit_opt$par[1],fit_opt$par[2],fit_opt$par[3],fit_opt$par[4])
  return(plot)
}
  

fit_hill_model_bc <- function(x_obs, y_obs, model, seed = 1) {
  hill_function_bc <- function(x, a, b,c, d) {
    a *(1 +(c-1)*(x^d / (b^d + x^d)))
  }
  
  neg_log_lik <- function(params, x, y) {
    a <- params[1]; b <- params[2]; c <- params[3];d <- params[4]; sigma <- params[5]
    y_hat <- hill_function_bc(x, a, b, c, d)
    -sum(dnorm(y, mean = y_hat, sd = sigma, log = TRUE))
  }
  
  runs <- map(1:5000, ~ {
    init <- c(
      a = runif(1, min(y_obs)/100, max(y_obs)*100), # baseline response
      # b = runif(1, 0.001, 100), # hill-ced model 
      b = runif(1, min(x_obs[x_obs > 0])/10, max(x_obs)*10), # for H5 model only,  inflection point (EC50 or BMD-like — dose at which half-max effect occurs)
      c = runif(1, 1, (max(y_obs)/min(y_obs))*5), # response scaling factor — asymptote at high dose, relative to baseline
      d = runif(1, 0.01, 10),  # slope or Hill coefficient — steepness of curve
      sigma = runif(1, 0.01, 0.1))
     optim(par = init, fn = neg_log_lik, x = x_obs, y = y_obs, method = "L-BFGS-B", 
          lower = c(min(y_obs)/100, min(x_obs[x_obs > 0])/10, 1, 0.01,1e-6), 
          upper = c(max(y_obs)*100, max(x_obs)*10, (max(y_obs)/min(y_obs))*5, 10, 1))
  })
  

  fit_opt <- runs[[which.min(purrr::map_dbl(runs, ~ .x$value))]]
  fit_opt$par
  
  stan_data <- list(
    N = length(x_obs),
    x = x_obs,
    y = y_obs,
    a_mean = max(fit_opt$par[1], 1e-6),
    a_sd = 0.1,
    b_mean = max(fit_opt$par[2], 1e-6),
    b_sd = 0.3, # 0.3
    c_mean = max(fit_opt$par[3], 1e-6),
    c_sd = 0.3, # 0.3
    d_mean = max(fit_opt$par[4], 1e-6),
    d_sd = 0.3 # 0.3
  )
  
  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 10000,
    iter_sampling = 10000,
    refresh = 0,
    show_messages = FALSE,
    adapt_delta = 0.999,
    max_treedepth = 80
  )
  
  return(fit)
}


plot_obs_pred <- function(x_obs,y_obs,a,b,c,d){
  x_lower <- if (min(x_obs) > 0) min(x_obs) / 10 else 1e-6
  x_upper <- max(x_obs) * 1
  x_seq <- exp(seq(log(x_lower), log(x_upper), length.out = 500))
  
  y_seq = hill_function_bc(x_seq,a,b,c,d)
  
  # Prepare data for the line
  line_df <- data.frame(x = x_seq, y = y_seq)
  
  # Plot
  ggplot() +
    geom_line(data = line_df, aes(x = x, y = y, color = "Median")) +
    geom_point(data = obs_df, aes(x = x, y = y, color = "Observed"), size = 2) +
    scale_color_manual(values = c("Median" = "blue", "Observed" = "red")) +
    labs(color = "Legend") +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = c(1, 1),
      legend.justification = c(0, 1),
      # legend.background = element_rect(fill = "white", color = "black"),
      legend.box.background = element_rect(color = "black")
    )
}


generate_pred_curves <- function(draws_df, x_seq, step = 100) {
  draw_indices <- seq(1, nrow(draws_df), by = step)
  
  pred_curves <- purrr::map_dfr(draw_indices, function(j) {
    a_j <- draws_df$a[j]
    b_j <- draws_df$b[j]
    c_j <- draws_df$c[j]
    d_j <- draws_df$d[j]
    
    y_pred <- a_j *(1 +(c_j-1)*(x_seq^d_j / (b_j^d_j + x_seq^d_j)))
    
    data.frame(x = x_seq, y = y_pred, draw = j)
  })
  
  ribbon_df <- pred_curves %>%
    dplyr::group_by(x) %>%
    dplyr::summarize(
      y_median = median(y),
      y_lower_95 = quantile(y, 0.025),
      y_upper_95 = quantile(y, 0.975),
      y_lower_50 = quantile(y, 0.25),
      y_upper_50 = quantile(y, 0.75)
    ) %>%
    dplyr::mutate(type_95 = "95% CI", type_50 = "50% CI")
  
  return(ribbon_df)
}


plot_hill_curve <- function(ribbon_df, obs_df, study_index, x_lower, x_upper,
                            ec = NULL,
                            ec_lower_50 = NULL, ec_upper_50 = NULL,
                            ec_lower_95 = NULL, ec_upper_95 = NULL,
                            y_lower = NULL,y_upper=NULL,yaxis_label = "Response") {

  p <- ggplot() +
    geom_ribbon(data = ribbon_df,
                aes(x = x, ymin = y_lower_95, ymax = y_upper_95, fill = type_95),
                alpha = 0.4) +
    geom_ribbon(data = ribbon_df,
                aes(x = x, ymin = y_lower_50, ymax = y_upper_50, fill = type_50),
                alpha = 0.6) +
    geom_line(data = ribbon_df,
              aes(x = x, y = y_median, color = "Median estimation"), size = 0.8) +
    geom_point(data = obs_df,
               aes(x = x, y = y, color = "Observed data"), size = 2) +
    scale_fill_manual(name = NULL,
                      values = c("95% CI" = "#C6DBEF", "50% CI" = "#6BAED6")) +
    scale_color_manual(name = NULL,
                       values = c("Observed data" = "black","Median estimation" = "#08519C"),
                       breaks = c("Observed data", "Median estimation")) +
    scale_x_log10(limits = c(x_lower, x_upper),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  oob = scales::oob_keep) +
    scale_y_continuous(limits = c(y_lower, y_upper)) +
    labs(
      #title = paste0("Study ", study_index, ": Predictive D-R Curve with EC"),
      x = expression(bolditalic("In vitro") ~ bold("concentration (μg/ml)")),
      y = yaxis_label,
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", size = 1,fill = NA),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 14,color = "black"),   # x-axis text size
      axis.title = element_text(size = 16, face = "bold",color = "black"),
      axis.ticks = element_line(color = "black"),
      #plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 12,color = "black"),
      legend.box.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
      legend.spacing.y = unit(-0.4, "cm")    # reduce vertical space between legends
    )
  
  # Add EC5 median vertical line

  # EC5 95% CI band
  if (!is.null(ec_lower_95) && !is.null(ec_upper_95)) {
    p <- p +
      annotate("rect", xmin = ec_lower_95, xmax = ec_upper_95,
               ymin = -Inf, ymax = Inf, alpha = 0.4, fill = "#F3A994")
  }
  #"#F9DBBE"
  # EC5 50% CI band
  if (!is.null(ec_lower_50) && !is.null(ec_upper_50)) {
    p <- p +
      annotate("rect", xmin = ec_lower_50, xmax = ec_upper_50,
               ymin = -Inf, ymax = Inf, alpha = 0.8, fill = "#F3A994")
  }
  
  # EC5 median vertical line
  if (!is.null(ec)) {
    y_pos <- if (!is.null(y_upper)) {
      (y_lower + y_upper) / 2
    } else {
      0.5 * (max(ribbon_df$y_upper_95, na.rm = TRUE) +
               min(ribbon_df$y_upper_95, na.rm = TRUE))
    }
    
    p <- p +
      geom_vline(xintercept = ec, linetype = "dashed", color = "#E9623C", size = 1) +
      annotate("text", x = ec_upper_50, y = y_pos,
               label = expression(bold(BMD[10])), vjust = -0.5, hjust = 0, color = "black",
               size = 5)
  }
  #"#C59368"

  return(p)
}

check_diagnostics <- function(fit, study_index, 
                              parameters = c("a", "b", "d","sigma", "ec5","ec10"),
                              rhat_threshold = 1.01, ess_threshold = 200) {
  summary_fit <- fit$summary()
  
  rhats <- summary_fit$rhat[summary_fit$variable %in% parameters]
  ess_bulk <- summary_fit$ess_bulk[summary_fit$variable %in% parameters]
  if (any(rhats > rhat_threshold, na.rm = TRUE)) {
    warning(sprintf("Study %d: Rhat > %.2f for some parameters", study_index, rhat_threshold))
  }
  
  if (any(ess_bulk < ess_threshold, na.rm = TRUE)) {
    warning(sprintf("Study %d: ESS < %d for some parameters", study_index, ess_threshold))
  }
}

save_traceplot <- function(fit, study_index, output_dir ){
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  traceplot <- bayesplot::mcmc_trace(
    posterior::as_draws_array(fit$draws()),
    pars = c("a", "b", "d", "sigma", "ec5","ec10")
  )
  
  ggsave(
    filename = sprintf("%s/study_%d_traceplot.pdf", output_dir, study_index),
    plot = traceplot,
    width = 8,
    height = 6
  )
}

model_bc <- cmdstan_model("with_stan/hill_bc_single.stan")


# --------------- 1. mouse biochemistry dataset unnormalized --------------------

load("/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/ls_sub_case_mouse.RData")
outputdir = "results_undissolved/mouse/plots/bmd_mouse_bc_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/stan_results/bmd_mouse_bc_stan"

ex_type = "biochemical"
col_response = 'Results'
outputfile = "bc_study_%d.rds"

ls_non_trend_human_bc = NaN


# ---------------- 2. human biochemistry dataset unnormalized -------------------

load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData")
outputdir = "results_undissolved/human/plots/bmd_human_bc_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results/bmd_human_bc_stan"


ex_type = "biochemical"
col_response = 'Results'
outputfile = "bc_study_%d.rds"

ls_non_trend_human_bc <- c(1, 4, 6, 10, 11, 12, 13, 19, 20, 21, 25, 26, 27, 29, 30, 35,
                           36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 
                           53,54,55,56,57,58, 76)

# ------------------ 3. human genototoxicity data  ------------------

load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData")
outputdir = "results_undissolved/human/plots/bmd_human_geno_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results/bmd_human_geno_stan"

ex_type = "genotoxicity"
col_response = 'Results'
outputfile = "geno_study_%d.rds"

ls_non_trend_human_bc = NaN

# ------------------4. human ldh data unnormalized -------

load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData")
outputdir = "results_undissolved/human/plots/bmd_human_ldh_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results/bmd_human_ldh_stan"


ex_type = "ldh"
col_response = 'Results'
outputfile = "ldh_study_%d.rds"

ls_non_trend_human_bc = NaN


# =========================== start of the code =================================

ls_bmd <- ls_sub_case[sapply(ls_sub_case, nrow) > 1]

ls_bmd_bc <- ls_bmd[grep(ex_type, names(ls_bmd), 
                           ignore.case = TRUE)]

unique_tested_assays <- unique(unlist(lapply(ls_bmd_bc, function(x) x$`Tested assay`)))
unique_tested_assays

# Clean and standardize assay names
library(stringr)

ls_bmd_bc <- lapply(ls_bmd_bc, function(df) {
  assays <- df$`Tested assay`
  
  # ---- Abbreviation: last (...) if present
  abbrevs <- str_extract(assays, "\\(([A-Za-z0-9\\-]+)\\)(?!.*\\([A-Za-z0-9\\-]+\\))")
  abbrevs <- str_remove_all(abbrevs, "[()]")
  
  # ---- Units: prefer last [ ... ], if not exist then last ( ... )
  units <- str_extract(assays, "\\[[^\\]]+\\](?!.*\\[[^\\]]+\\])")
  units <- str_remove_all(units, "[\\[\\]]")
  
  missing_units <- is.na(units) | units == ""
  units[missing_units] <- str_extract(assays[missing_units], "\\([^()]+\\)(?!.*\\([^()]+\\))")
  units <- str_remove_all(units, "[()]")
  
  # ---- Remove leading category prefixes
  clean_name <- str_remove(assays, "^(Biochemical:|Cytokine:)\\s*")
  clean_name <- str_trim(str_remove(clean_name, "\\[.*"))  # drop [..] from description
  
  # ---- Final label ----
  final <- ifelse(
    is.na(abbrevs) | abbrevs == "",
    paste0(clean_name, " [", units, "]"),   # descriptive name + unit
    paste0(abbrevs, " [", units, "]")       # abbreviation + unit
  )
  
  df$`Tested assay` <- final
  df
})

ls_bmd_bc <- lapply(ls_bmd_bc, function(df) {
  assays <- df$`Tested assay`
  
  # Special replacement for ROS entries
  assays <- gsub("ROS \\(Oxidative damage\\) \\[DCFH-DA assay\\]",
                 "DCFH-DA assay [% of ROS production]",
                 assays)
  
  df$`Tested assay` <- assays
  df
})


# Process y-values for each element in ls_bmd_bc based on assay type (normalization for different assays)
for (i in seq_along(ls_bmd_bc)) {
  assay <- unique(ls_bmd_bc[[i]]$`Tested assay`)
  y <- ls_bmd_bc[[i]]$Results
  control_value <- y[1]
  
  # Initialize Processed column
  ls_bmd_bc[[i]]$Processed <- y
  
  # Case 1: GSH given as % of control → compute depletion (e.g., GSH depletion)
  if (grepl("GSH", assay, ignore.case = TRUE) && grepl("% of control", assay)) {
    ls_bmd_bc[[i]]$Processed <- pmax(100 - y, 1e-6)
    
    # Case 2: MDA or similar markers as % of control → increase above 100%
  } else if (grepl("MDA|CAT", assay, ignore.case = TRUE) && grepl("% of control", assay)) {
    ls_bmd_bc[[i]]$Processed <- pmax(y - 100, 1e-6)
    
    # Default case: Normalize to control value
  } else if (grepl("GSH|protein assay|GPx|thrr|G6PDH|GRD|SOD", assay, ignore.case = TRUE) && !grepl("% of control", assay)) {
    ls_bmd_bc[[i]]$Processed <- pmax(y / control_value * 100, 1e-6)
    ls_bmd_bc[[i]]$Processed <- pmax(100 - ls_bmd_bc[[i]]$Processed, 1e-6)
  } 
  
  else {
    ls_bmd_bc[[i]]$Processed <- pmax((y-control_value) / control_value * 100, 1e-6)
  }
  
}

# ----------- start of the dose response curve fitting --------------

results <- data.frame(study = integer(), a = numeric(), b = numeric(), d = numeric(),
                      sigma = numeric(), ec5_median = numeric(), ec5_q5 = numeric(), 
                      ec5_q25 = numeric(), ec5_q75 = numeric(), ec5_q95 = numeric(),
                      ec10_median = numeric(), ec10_q5 = numeric(), ec10_q25 = numeric(),
                      ec10_q75 = numeric(), ec10_q95 = numeric(),r2 = numeric())
ec5_draws_all <- list()
ec10_draws_all <- list()
param_draws_all <- list()
r_squared_df <- list()

for (i in seq_along(ls_bmd_bc)) {

  percent <- round(100 * i / length(ls_bmd_bc))
  cat(sprintf("\rProgress: %3d%%", percent))
  flush.console()
  
  data_i <- ls_bmd_bc[[i]]
  x_obs <- data_i$`In vitro concentration`
  y_obs <- data_i[,col_response]  # Use normalized results or Use raw values
  
  obs_df <- data.frame(x = x_obs, y = y_obs)
  obs_df
  y_axis_title <- data_i$`Tested assay` %>%
    unique() %>%
    str_remove("^Biochemical:\\s*")%>%
    str_remove("^Genotoxicity:\\s*")%>%
    str_remove("\\(DNA damage\\) ") %>%   # remove "(DNA damage)"
    str_replace(" to negative control", "")
  
  if (nrow(obs_df) == 0) {next}
  if (all(y_obs == 1e-6)){   
    p = ggplot() +
      geom_point(data = data_i, aes(x = `In vitro concentration`, y = Results)
                 ,color="black", size = 2) +
      theme_minimal(base_size = 12) +
      theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5),
      )
    ggsave(
      file.path(outputdir, sprintf("opposite_trend/study_%d_curve.pdf", i)),
      plot = p, width = 5, height = 4) 
    next}
  if (i %in% ls_non_trend_human_bc){
    p = ggplot() +
    geom_point(data = data_i, aes(x = `In vitro concentration`, y = Results)
               ,color="black", size = 2) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
    )
  ggsave(
    file.path(outputdir, sprintf("no_trend/study_%d_curve.pdf", i)),
    plot = p, width = 5, height = 4) 
  next}
  
  obs_df <- obs_df[!(obs_df$y == 1e-6 & seq_len(nrow(obs_df)) != 1), ] # removes all rows where `y == 1e-6`**, **except the first row
  
  if (nrow(obs_df) == 2 & min(obs_df)==0) {
    p = plot_fit_initial(x_obs,y_obs)
    ggsave(
      file.path(outputdir, sprintf("single/study_%d_curve.pdf", i)),
      plot = p, width = 8, height = 6) 
    next} # skip if only one data points are available 
  
  # fit the model and then save the results 
  fit <- fit_hill_model_bc(obs_df$x, obs_df$y, model_bc)  # Fit model
  dir.create(dirname(file.path(stanoutdir, sprintf(outputfile, i))), recursive = TRUE, showWarnings = FALSE)
  
  fit$save_object(file = file.path(stanoutdir, sprintf(outputfile, i)))   # Save a single object

  # Summarize results

  draws_df <- posterior::as_draws_df(fit$draws())
  
  summary_pars <- summarise_draws(
    draws_df[, c("a", "b", "c","d", "sigma", "ec5", "ec10")],   # your posterior draws
    mean, median, sd, mad,
    ~quantile2(.x, probs = c(0.05, 0.25, 0.75, 0.95)),  # add q25 here
    rhat, ess_bulk, ess_tail
  )
  
  # Predictive curves


  median_ec5 <- median(draws_df$ec5)
  ec5_ci_50 <- quantile(draws_df$ec5, probs = c(0.25, 0.75))
  ec5_ci_95 <- quantile(draws_df$ec5, probs = c(0.025, 0.975))
  
  median_ec10 <- median(draws_df$ec10)
  ec10_ci_50 <- quantile(draws_df$ec10, probs = c(0.25, 0.75))
  ec10_ci_95 <- quantile(draws_df$ec10, probs = c(0.025, 0.975))
  


  x_lower <- if (min(x_obs) > 0) min(x_obs) / 10 else 1e-3
  x_lower <- if (x_lower > ec10_ci_95[1]) ec10_ci_95[1]/10 else x_lower
  x_upper <- max(x_obs) * 5 # for the cases that x axis should be zoomed 
  x_seq <- exp(seq(log(x_lower), log(x_upper), length.out = 500))
  ribbon_df <- generate_pred_curves(draws_df, x_seq)
  
  p <- plot_hill_curve(ribbon_df, obs_df, i, x_lower, x_upper, 
                       median_ec10,ec10_ci_50[1],ec10_ci_50[2],
                       ec10_ci_95[1],ec10_ci_95[2],
                       yaxis_label=capitalize_if_needed(y_axis_title))
  p
  dir.create(dirname(file.path(outputdir, sprintf("study_%d_curve.pdf", i))), recursive = TRUE, showWarnings = FALSE)
  
  if (abs(log10(median_ec10) - log10(ec10_ci_95[2])) > 3) {
    # Save plot
    ggsave(
      file.path(outputdir, sprintf("uncertain/study_%d_curve.pdf", i)),
      plot = p, width = 8, height = 5,device = cairo_pdf)}else{  # Save plot
    ggsave(
         file.path(outputdir, sprintf("study_%d_curve.pdf", i)),
         plot = p, width =8, height = 5,device = cairo_pdf)}
  
  # After sampling and summarizing
  check_diagnostics(fit, study_index = i)
  #save_traceplot(fit, study_index = i,outputdir)
  
  ec5_draws_all[[i]] <- data.frame(
    study = paste0("Study_", i),
    ec5 = draws_df$ec5
  )
  
  ec10_draws_all[[i]] <- data.frame(
    study = paste0("Study_", i),
    ec10 = draws_df$ec10
  )
  
  param_draws_all[[i]] <- subset_draws(draws_df, variable = c("a", "b", "d")) %>%
    as_draws_df() %>%
    as.data.frame() %>%
    dplyr::select(a, b, d) %>%   # explicitly drop .chain, etc.
    dplyr::mutate(study = paste0("Study_", i))
  
  
  # calculate the R square and the observation vs prediction
  # For each x in obs_df, find closest x in ribbon_df
  matched_preds <- sapply(obs_df$x, function(x_obs) {
    ribbon_df$y_median[which.min(abs(ribbon_df$x - x_obs))]
  })
  
  # Create aligned data frame
  aligned_df <- data.frame(
    x = obs_df$x,
    y_obs = obs_df$y,
    y_pred = matched_preds
  )
  
  # Compute R-squared
  ss_res <- sum((aligned_df$y_obs - aligned_df$y_pred)^2, na.rm = TRUE)
  ss_tot <- sum((aligned_df$y_obs - mean(aligned_df$y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- 1 - ss_res / ss_tot
  
  aligned_df$r_squared = r_squared
  
  r_squared_df[[i]] = aligned_df
  rm(fit)
  
  results <- dplyr::bind_rows(results, data.frame(
    study = i,
    a = summary_pars$mean[summary_pars$variable == "a"],
    b = summary_pars$mean[summary_pars$variable == "b"],
    c = summary_pars$mean[summary_pars$variable == "c"],
    d = summary_pars$mean[summary_pars$variable == "d"],
    sigma = summary_pars$mean[summary_pars$variable == "sigma"],
    ec5_median = summary_pars$median[summary_pars$variable == "ec5"],
    ec5_q5 = summary_pars$q5[summary_pars$variable == "ec5"],
    ec5_q25 = summary_pars$q25[summary_pars$variable == "ec5"],
    ec5_q75 = summary_pars$q75[summary_pars$variable == "ec5"],
    ec5_q95 = summary_pars$q95[summary_pars$variable == "ec5"],
    ec10_median = summary_pars$median[summary_pars$variable == "ec10"],
    ec10_q5 = summary_pars$q5[summary_pars$variable == "ec10"],
    ec10_q25 = summary_pars$q25[summary_pars$variable == "ec10"],
    ec10_q75 = summary_pars$q75[summary_pars$variable == "ec10"],
    ec10_q95 = summary_pars$q95[summary_pars$variable == "ec10"],
    r2 = r_squared
  ))
  
}

# ============================ save the results =================================

# --------------- 2. mouse biochemistry dataset unnormalized --------------------
save(param_draws_all, file = file.path(stanoutdir,"bc_mouse_param_draws_all.RData"))
save(ec5_draws_all, file =file.path(stanoutdir,"bc_mouse_ec5_draws_all.RData"))
save(ec10_draws_all, file =file.path(stanoutdir,"bc_mouse_ec10_draws_all.RData"))
save(results, file= file.path(stanoutdir,"bc_mouse_results_summary.RData"))
save(r_squared_df, file = file.path(stanoutdir,"bc_mouse_r_square.RData"))


# ---------------- 4. human biochemistry dataset unnormalized -------------------
save(param_draws_all, file = file.path(stanoutdir,"bc_human_param_draws_all.RData"))
save(ec5_draws_all, file =file.path(stanoutdir,"bc_human_ec5_draws_all.RData"))
save(ec10_draws_all, file =file.path(stanoutdir,"bc_human_ec10_draws_all.RData"))
save(results, file= file.path(stanoutdir,"bc_human_results_summary.RData"))
save(r_squared_df, file = file.path(stanoutdir,"bc_human_r_square.RData"))



# ------------------ 6. human genototoxicity data unnormalized ------------------
save(param_draws_all, file = file.path(stanoutdir,"geno_human_param_draws_all.RData"))
save(ec5_draws_all, file =file.path(stanoutdir,"geno_human_ec5_draws_all.RData"))
save(ec10_draws_all, file =file.path(stanoutdir,"geno_human_ec10_draws_all.RData"))
save(results, file= file.path(stanoutdir,"geno_human_results_summary.RData"))
save(r_squared_df, file = file.path(stanoutdir,"geno_human_r_square.RData"))

# ------------------ 7. human ldh data unnormalized ------------------
save(param_draws_all, file = file.path(stanoutdir,"ldh_human_param_draws_all.RData"))
save(ec5_draws_all, file =file.path(stanoutdir,"ldh_human_ec5_draws_all.RData"))
save(ec10_draws_all, file =file.path(stanoutdir,"ldh_human_ec10_draws_all.RData"))
save(results, file= file.path(stanoutdir,"ldh_human_results_summary.RData"))
save(r_squared_df, file = file.path(stanoutdir,"ldh_human_r_square.RData"))


# ======= After sampling and summarizing, checking if everyone converges ========
for (i in seq_along(ls_bmd_bc)) {
#for (i in ls_non_converge){
  percent <- round(100 * i / length(ls_bmd_bc))
  cat(sprintf("\rProgress: %3d%%", percent))
  flush.console()

  fit <- tryCatch({
    readRDS(file.path(stanoutdir, sprintf(outputfile, i)))
  }, error = function(e) {
    message(sprintf("Skipping study %d due to error: %s", i, e$message))
    return(NULL)  # or NA, or next
  })
  
  #check_diagnostics(summary_fit, study_index = i)
  if (!is.null(fit)) {
    check_diagnostics(fit, study_index = i)
  }
}

# for human biochemistry results 
# the non calculated list: usually it is due to the results are not significant 
# different from the benchmark value, then after normalization, everything is smaller than 0
ls_non_converge = c(6,63,64,73,85,128,141,145,155,171,174,176,185,189,196,197) 
ls_non_calculate = c(25,30,36,39,42,54,55,57,110,115,125,136,148) 

ls_non_converge = c(2,7,14,155,196,200)
ls_onepoint = c(11,24,25,30,36,39,42,44,54,55,57,59,65,69,70,71,86,89,96,101,104,
                110,113,115,122,123,125,136,137,142,143,148,158,163,164,201)



ec5_draws_df <- dplyr::bind_rows(ec5_draws_all)

# Reorder by EC5 median (optional)
study_order <- ec5_draws_df %>%
  group_by(study) %>%
  summarize(median_ec5 = median(ec5)) %>%
  arrange(median_ec5) %>%
  pull(study)

ec5_draws_df$study <- factor(ec5_draws_df$study, levels = study_order)

# Plot
ggplot(ec5_draws_df, aes(x = study, y = ec5)) +
  geom_violin(fill = "#DDEEFF", color = "black", scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_minimal(base_size = 12) +
  scale_y_log10()+
  labs(title = "Posterior EC5 distributions by study",
       x = "Study",
       y = "EC5 (Effective Concentration)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/mouse/plots/ec5_distributions_violin.pdf", width = 10, height = 6)


