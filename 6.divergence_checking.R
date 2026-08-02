"results_undissolved/mouse/stan_results/bmd_mouse_bc_stan/bc_mouse_diagnostics.RData"
"results_undissolved/mouse/stan_results/bmd_mouse_cyto_stan/cyto_mouse_diagnostics.RData"
"results_undissolved/human/stan_results/bmd_human_bc_stan/bc_human_diagnostics.RData"
"results_undissolved/human/stan_results/bmd_human_cyto_stan/cyto_human_diagnostics.RData"
"results_undissolved/human/stan_results/bmd_human_geno_stan/geno_human_diagnostics.RData"
"results_undissolved/human/stan_results/bmd_human_ldh_stan/ldh_human_diagnostics.RData"

library(dplyr)

files <- c(
  "results_undissolved/mouse/stan_results/bmd_mouse_bc_stan/bc_mouse_diagnostics.RData",
  "results_undissolved/mouse/stan_results/bmd_mouse_cyto_stan/cyto_mouse_diagnostics.RData",
  "results_undissolved/human/stan_results/bmd_human_bc_stan/bc_human_diagnostics.RData",
  "results_undissolved/human/stan_results/bmd_human_cyto_stan/cyto_human_diagnostics.RData",
  "results_undissolved/human/stan_results/bmd_human_geno_stan/geno_human_diagnostics.RData",
  "results_undissolved/human/stan_results/bmd_human_ldh_stan/ldh_human_diagnostics.RData"
)

source_name <- c(
  "Matabolic cell stress - mouse",
  "Cytotoxicity - mouse",
  "Matabolic cell stress - human",
  "Cytotoxicity - human",
  "Genotoxicity - human",
  "LDH - human"
)

stopifnot(length(files) == length(source_name))   # sanity check

load_diagnostics <- function(file, label) {
  e <- new.env()
  obj_names <- load(file, envir = e)
  df <- as.data.frame(e[[obj_names[1]]])
  df$source <- label
  names(df)[names(df) == "study"] <- "CASE_ID"
  df
}

merged_diagnostics <- dplyr::bind_rows(Map(load_diagnostics, files, source_name))

# Save the merged result
readr::write_csv(merged_diagnostics, "results_undissolved/merged_diagnostics.csv")


