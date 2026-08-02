# ==============================================================================
#
# TITLE:   Dose-response curve fitting for cytotoxicity experiments
#
# PURPOSE: Takes the cleaned, per-experiment in vitro data and fits dose-response
#          curves for the cytotoxicity endpoint class, covering cell viability
#          assays (MTT, Alamar Blue, EZ-Cytox, MTS, WST-1, Trypan Blue, NRU,
#          ATP, CCK8) and the membrane integrity assay (LDH).
#          These endpoints have a biologically interpretable response range
#          (100% viability -> 0%), so a three-parameter Hill function with a
#          fixed upper asymptote is used:
#
#              y = a + (1 - a) * x^d / (b^d + x^d)
#
#          where x = nanoparticle dose, y = normalized viability, a = baseline
#          viability, b = half-maximal dose, d = slope; maximum fixed at 1.
#          LDH is handled on the inverse (% cytotoxicity) scale.
#          Fitting is Bayesian (MCMC in R); convergence checked via R-hat < 1.2.
#          BMD10 (10% change from the negative control) is extracted per
#          experiment together with its credible interval.
#
# INPUT:   <cleaned per-experiment data set from 01_read_clean_invitro_data.R>
# OUTPUT:  <fitted model objects / BMD10 table + dose-response curve plots>
#
# AUTHOR:  Jimeng Wu          CREATED: 2026-04        
# ==============================================================================



library(posterior)
library(purrr)
library(cmdstanr)
library(ggplot2)
library(dplyr)

fit_hill_model <- function(x_obs, y_obs, model, seed = 234) {
  hill_function <- function(x, a, b, d) {
    a + (1 - a) * (x^d / (b^d + x^d))
  }
  
  neg_log_lik <- function(params, x, y) {
    a <- params[1]; b <- params[2]; d <- params[3]; sigma <- params[4]
    y_hat <- hill_function(x, a, b, d)
    -sum(dnorm(y, mean = y_hat, sd = sigma, log = TRUE))
  }
  
  runs <-map(1:2000, ~ {
    init <- c(
      a = runif(1, 0.001, 0.01), 
      b = runif(1, 0.001, 10000), 
      d = runif(1, 0.01, 10), 
      sigma = runif(1, 0.01, 0.1))
    optim(par = init, fn = neg_log_lik, x = x_obs, y = y_obs, method = "L-BFGS-B", 
          lower = c(1e-06, 1e-3, 0.1,1e-6), upper = c(0.01, 1e+06, 10, 1))
  })
  
  
 
  
  fit_opt <- runs[[which.min(purrr::map_dbl(runs, ~ .x$value))]]
  
  stan_data <- list(
    N = length(x_obs),
    x = x_obs,
    y = y_obs,
    a_mean = max(fit_opt$par[1], 1e-6),
    a_sd = 0.1,
    b_mean = max(fit_opt$par[2], 1e-6),
    b_sd = 0.3,
    d_mean = max(fit_opt$par[3], 1e-6),
    d_sd = 0.3
  )
  
  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 10000,
    iter_sampling = 10000,
    refresh = 500,
    adapt_delta = 0.995,
    max_treedepth = 15
  )
  
  return(fit)
}


generate_pred_curves <- function(draws_df, x_seq, step = 100) {
  draw_indices <- seq(1, nrow(draws_df), by = step)
  
  pred_curves <- purrr::map_dfr(draw_indices, function(j) {
    a_j <- draws_df$a[j]
    b_j <- draws_df$b[j]
    d_j <- draws_df$d[j]
    
    y_pred <- a_j + (1 - a_j) * (x_seq^d_j / (b_j^d_j + x_seq^d_j))
    
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


model <- cmdstan_model("with_stan/hill_single.stan")

# mouse cyto toxicity dataset 
load("/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/ls_sub_case_mouse.RData")
outputdir = "results_undissolved/mouse/plots/bmd_mouse_cyto_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/stan_results/bmd_mouse_cyto_stan"
plot_title_begin="Cytotoxicity"

# human cyto toxicity dataset
load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData")
outputdir = "results_undissolved/human/plots/bmd_human_cyto_stan"
stanoutdir = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results/bmd_human_cyto_stan"


#====================== start of the code =========================
ls_bmd <- ls_sub_case[sapply(ls_sub_case, nrow) > 1]
ls_bmd_cyto <- ls_bmd[grep("% cytotoxicity|viability % of control", names(ls_bmd), 
                           ignore.case = TRUE)]

unique_tested_assays <- unique(unlist(lapply(ls_bmd_cyto, function(x) x$`Tested assay`)))
unique_tested_assays


results <- data.frame(study = integer(), a = numeric(), b = numeric(), d = numeric(),
                      sigma = numeric(), ec5_median = numeric(), ec5_q5 = numeric(), 
                      ec5_q25 = numeric(), ec5_q75 = numeric(), ec5_q95 = numeric(),
                      ec10_median = numeric(), ec10_q5 = numeric(), ec10_q25 = numeric(),
                      ec10_q75 = numeric(), ec10_q95 = numeric(),r2 = numeric())

ec5_draws_all <- list()
ec10_draws_all <- list()
param_draws_all <- list()
r_squared_df <- list()
all_diagnostics <- vector("list", length(ls_bmd_cyto))

for (i in seq_along(ls_bmd_cyto)) {
  percent <- round(100 * i / length(ls_bmd_cyto))
  cat(sprintf("\rProgress: %3d%%", percent))
  flush.console()
  
  data_i <- ls_bmd_cyto[[i]]
  x_obs <- data_i$`In vitro concentration`
  y_obs <- data_i$Results / 100
  obs_df <- data.frame(x = x_obs, y = y_obs)
  obs_df
  
  #  1. fit the model and then save the results 
  #fit <- fit_hill_model(x_obs, y_obs, model) # Fit model
  #fit$save_object(file = file.path(stanoutdir, sprintf("cyto_study_%d.rds", i)))   # Save a single object
  # 2. read the saved model fitting results
  fit <- readRDS(file.path(stanoutdir, sprintf("cyto_study_%d.rds", i)))
  
  
  # Summarize results
  #summary_pars <- fit$summary(c("a", "b", "d", "sigma", "ec5","ec10"))
  draws_df <- posterior::as_draws_df(fit$draws())
  
  summary_pars <- summarise_draws(
    draws_df[, c("a", "b", "d", "sigma", "ec5", "ec10")],   # your posterior draws
    mean, median, sd, mad,
    ~quantile2(.x, probs = c(0.05, 0.25, 0.75, 0.95)),  # add q25 here
    rhat, ess_bulk, ess_tail
  )
  

  

  median_ec10 <- median(draws_df$ec10)
  ec10_ci_50 <- quantile(draws_df$ec10, probs = c(0.25, 0.75))
  ec10_ci_95 <- quantile(draws_df$ec10, probs = c(0.025, 0.975))
  
  # Predictive curves
  x_lower <- if (min(x_obs) > 0) min(x_obs) / 10 else 1e-3
  x_lower <- if (x_lower > ec10_ci_95[1]) ec10_ci_95[1]/10 else x_lower
  x_upper <- max(x_obs) * 100
  x_upper <- if (x_upper < ec10_ci_95[2]) ec10_ci_95[2]*10 else x_upper
  
  x_seq <- exp(seq(log(if (min(x_obs) > 0) x_lower else 1e-6), log(x_upper), length.out = 500))
  
  ribbon_df <- generate_pred_curves(draws_df, x_seq)
  p <- plot_hill_curve(ribbon_df, ls_bmd_cyto, obs_df, i, x_lower, x_upper, 
                       median_ec10,ec10_ci_50[1],ec10_ci_50[2],
                       ec10_ci_95[1],ec10_ci_95[2],0,1,"Cytotoxicity (%)")
  
  # Save plot
  ggsave(
        file.path(outputdir, sprintf("study_%d_curve.pdf", i)),
       plot = p, width = 8, height = 5,device = cairo_pdf)
  
  # After sampling and summarizing
  all_diagnostics[[i]] <- check_diagnostics(fit, study_index = i)
  
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
  
  results <- dplyr::bind_rows(results, data.frame(
    study = i,
    a = summary_pars$mean[summary_pars$variable == "a"],
    b = summary_pars$mean[summary_pars$variable == "b"],
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
diagnostics_df <- do.call(rbind, all_diagnostics)

# for mouse
#save(param_draws_all, file = file.path(stanoutdir,"cyto_mouse_param_draws_all.RData"))
#save(ec5_draws_all, file = file.path(stanoutdir,"cyto_mouse_ec5_draws_all.RData"))
#save(ec10_draws_all, file = file.path(stanoutdir,"cyto_mouse_ec10_draws_all.RData"))
#save(results, file = file.path(stanoutdir,"cyto_mouse_results_summary.RData"))
#save(r_squared_df, file = file.path(stanoutdir,"cyto_mouse_r_square.RData"))

save(diagnostics_df, file = file.path(stanoutdir,"cyto_mouse_diagnostics.RData"))

# for human
#save(param_draws_all, file = file.path(stanoutdir,"cyto_human_param_draws_all.RData"))
#save(ec5_draws_all, file = file.path(stanoutdir,"cyto_human_ec5_draws_all.RData"))
#save(ec10_draws_all, file = file.path(stanoutdir,"cyto_human_ec10_draws_all.RData"))
#save(results, file = file.path(stanoutdir,"cyto_human_results_summary.RData"))
#save(r_squared_df, file = file.path(stanoutdir,"cyto_human_r_square.RData"))

save(diagnostics_df, file = file.path(stanoutdir,"cyto_human_diagnostics.RData"))


# After sampling and summarizing
for (i in seq_along(ls_bmd_cyto)) {
  percent <- round(100 * i / length(ls_bmd_cyto))
  cat(sprintf("\rProgress: %3d%%", percent))
  flush.console()
  fit <- readRDS(file.path(stanoutdir, sprintf("cyto_study_%d.rds", i)))
  
  check_diagnostics(fit, study_index = i)
}




