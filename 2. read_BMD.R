library(dplyr)
library(stringr)
library(scales)
library(forcats)
source("with_stan/helper_function.R")


stanoutdir_mouse = "/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/stan_results"
stanoutdir_human = "/Users/wuji/work/code/codo_v2/R/results_undissolved/human/stan_results"

# for human cytotoxicity
cyto_human_param_draws_all = get(load(file = file.path(stanoutdir_human,"bmd_human_cyto_stan/cyto_human_param_draws_all.RData")))
cyto_human_ec5_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_cyto_stan/cyto_human_ec5_draws_all.RData")))
cyto_human_ec10_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_cyto_stan/cyto_human_ec10_draws_all.RData")))
cyto_human_results = get(load(file= file.path(stanoutdir_human,"bmd_human_cyto_stan/cyto_human_results_summary.RData")))
cyto_human_r_sqaure = get(load(file = file.path(stanoutdir_human,"bmd_human_cyto_stan/cyto_human_r_square.RData")))

# for mouse cytotoxicity
cyto_mouse_param_draws_all = get(load(file = file.path(stanoutdir_mouse,"bmd_mouse_cyto_stan/cyto_mouse_param_draws_all.RData")))
cyto_mouse_ec5_draws_all = get(load(file =file.path(stanoutdir_mouse,"bmd_mouse_cyto_stan/cyto_mouse_ec5_draws_all.RData")))
cyto_mouse_ec10_draws_all = get(load(file =file.path(stanoutdir_mouse,"bmd_mouse_cyto_stan/cyto_mouse_ec10_draws_all.RData")))
cyto_mouse_results = get(load(file= file.path(stanoutdir_mouse,"bmd_mouse_cyto_stan/cyto_mouse_results_summary.RData")))
cyto_mouse_r_sqaure = get(load(file = file.path(stanoutdir_mouse,"bmd_mouse_cyto_stan/cyto_mouse_r_square.RData")))

# for human biochemistry
bc_human_param_draws_all = get(load(file = file.path(stanoutdir_human,"bmd_human_bc_stan/bc_human_param_draws_all.RData")))
bc_human_ec5_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_bc_stan/bc_human_ec5_draws_all.RData")))
bc_human_ec10_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_bc_stan/bc_human_ec10_draws_all.RData")))
bc_human_results = get(load(file= file.path(stanoutdir_human,"bmd_human_bc_stan/bc_human_results_summary.RData")))
bc_human_r_sqaure = get(load(file = file.path(stanoutdir_human,"bmd_human_bc_stan/bc_human_r_square.RData")))

# for mouse biochemistry
bc_mouse_param_draws_all = get(load(file = file.path(stanoutdir_mouse,"bmd_mouse_bc_stan/bc_mouse_param_draws_all.RData")))
bc_mouse_ec5_draws_all = get(load(file =file.path(stanoutdir_mouse,"bmd_mouse_bc_stan/bc_mouse_ec5_draws_all.RData")))
bc_mouse_ec10_draws_all = get(load(file =file.path(stanoutdir_mouse,"bmd_mouse_bc_stan/bc_mouse_ec10_draws_all.RData")))
bc_mouse_results = get(load(file= file.path(stanoutdir_mouse,"bmd_mouse_bc_stan/bc_mouse_results_summary.RData")))
bc_mouse_r_sqaure = get(load(file = file.path(stanoutdir_mouse,"bmd_mouse_bc_stan/bc_mouse_r_square.RData")))

# for human genotoxicity
geno_human_param_draws_all = get(load(file = file.path(stanoutdir_human,"bmd_human_geno_stan/geno_human_param_draws_all.RData")))
geno_human_ec5_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_geno_stan/geno_human_ec5_draws_all.RData")))
geno_human_ec10_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_geno_stan/geno_human_ec10_draws_all.RData")))
geno_human_results = get(load(file= file.path(stanoutdir_human,"bmd_human_geno_stan/geno_human_results_summary.RData")))
geno_human_r_sqaure = get(load(file = file.path(stanoutdir_human,"bmd_human_geno_stan/geno_human_r_square.RData")))


# for human ldh
ldh_human_param_draws_all = get(load(file = file.path(stanoutdir_human,"bmd_human_ldh_stan/ldh_human_param_draws_all.RData")))
ldh_human_ec5_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_ldh_stan/ldh_human_ec5_draws_all.RData")))
ldh_human_ec10_draws_all = get(load(file =file.path(stanoutdir_human,"bmd_human_ldh_stan/ldh_human_ec10_draws_all.RData")))
ldh_human_results = get(load(file= file.path(stanoutdir_human,"bmd_human_ldh_stan/ldh_human_results_summary.RData")))
ldh_human_r_sqaure = get(load(file = file.path(stanoutdir_human,"bmd_human_ldh_stan/ldh_human_r_square.RData")))


# for mouse cytokine
cytokine_mouse_ec5_draws_all = get(load(file=file.path(stanoutdir_mouse, "bmd_mouse_cytokine_stan/cytokine_mouse_ec5_draws_all.RData")))

# for human cytokine
cytokine_human_ec5_draws_all = get(load(file=file.path(stanoutdir_human, "bmd_human_cytokine_stan/cytokine_human_ec5_draws_all.RData")))



# ======================2. get the dataset inforamtion ===================

ls_sub_case_human = get(load("/Users/wuji/work/code/codo_v2/R/results/human/ls_sub_case_human.RData"))

ls_sub_case_human = get(load("/Users/wuji/work/code/codo_v2/R/results_undissolved/human/ls_sub_case_human.RData"))

ls_bmd_human <- ls_sub_case_human[sapply(ls_sub_case_human, nrow) > 1]

ls_bmd_cyto_human <- ls_bmd_human[grep("% cytotoxicity|viability % of control", names(ls_bmd_human), 
                           ignore.case = TRUE)]

ls_bmd_bc_human <- ls_bmd_human[grep("biochemical", names(ls_bmd_human), 
                         ignore.case = TRUE)]

ls_bmd_cytokine_human <- ls_bmd_human[grep("cytokine", names(ls_bmd_human), 
                               ignore.case = TRUE)]

ls_bmd_geno_human <- ls_bmd_human[grep("genotoxicity", names(ls_bmd_human), 
                         ignore.case = TRUE)]

ls_bmd_ldh_human <- ls_bmd_human[grep("ldh", names(ls_bmd_human), 
                                       ignore.case = TRUE)]

length(ls_sub_case_human)
length(ls_bmd_human)
length(ls_bmd_cyto_human) +length(ls_bmd_bc_human) + length(ls_bmd_cytokine_human)+length(ls_bmd_geno_human) + length(ls_bmd_ldh_human)


# mouse
ls_sub_case_mouse = get(load("/Users/wuji/work/code/codo_v2/R/results/mouse/ls_sub_case_mouse.RData"))
ls_sub_case_mouse = get(load("/Users/wuji/work/code/codo_v2/R/results_undissolved/mouse/ls_sub_case_mouse.RData"))
length(ls_sub_case_mouse)
ls_bmd_mouse <- ls_sub_case_mouse[sapply(ls_sub_case_mouse, nrow) > 1]
length(ls_bmd_mouse)
ls_bmd_cyto_mouse <- ls_bmd_mouse[grep("% cytotoxicity|viability % of control", names(ls_bmd_mouse), 
                                       ignore.case = TRUE)]

ls_bmd_bc_mouse <- ls_bmd_mouse[grep("biochemical", names(ls_bmd_mouse), 
                                     ignore.case = TRUE)]

ls_bmd_cytokine_mouse <- ls_bmd_mouse[grep("cytokine", names(ls_bmd_mouse), 
                                           ignore.case = TRUE)]
length(ls_bmd_cyto_mouse) +length(ls_bmd_bc_mouse) + length(ls_bmd_cytokine_mouse)

# ----- in paper we said 144 for human, but quantum dot was removed here
# ----------------------- 2.1 check if it covers all -----------------------

# Get all names from the three sublists
all_names_combined_human <- c(
  names(ls_bmd_cytokine_human),
  names(ls_bmd_bc_human),
  names(ls_bmd_cyto_human),
  names(ls_bmd_geno_human),
  names(ls_bmd_ldh_human)
)

# Get the unique names in the union of the three sublists

# Get the names in the main list
main_names_human <- names(ls_bmd_human)

# Check if all names in ls_bmd_human are covered
all_covered_human <- all(names(ls_bmd_human) %in% unique(all_names_combined_human))

# Show uncovered names if any
uncovered_names_human <- setdiff(main_names_human, unique(all_names_combined_human))

# Output
cat("All items for human covered:", all_covered_human, "\n")
if (!all_covered_human) {
  cat("Uncovered items:\n")
  print(uncovered_names_human)
}


# or mouse
all_names_combined_mouse <- c(
  names(ls_bmd_cytokine_mouse),
  names(ls_bmd_bc_mouse),
  names(ls_bmd_cyto_mouse)
)
main_names_mouse <- names(ls_bmd_mouse)
# Check if all names in ls_bmd_human are covered
all_covered_mouse <- all(names(ls_bmd_mouse) %in% unique(all_names_combined_mouse))

# Show uncovered names if any
uncovered_names_mouse<- setdiff(main_names_mouse, unique(all_names_combined_mouse))

# Output
cat("All items for mouse covered:", all_covered_mouse, "\n")
if (!all_covered_mouse) {
  cat("Uncovered items:\n")
  print(uncovered_names_mouse)
}


# =====================3. match the experiment information to the ec results ===========


cyto_mouse_combined <- generate_combined_results(ls_bmd_cyto_mouse, cyto_mouse_results)
cyto_human_combined <- generate_combined_results(ls_bmd_cyto_human, cyto_human_results)
bc_mouse_combined   <- generate_combined_results(ls_bmd_bc_mouse, bc_mouse_results)
bc_human_combined   <- generate_combined_results(ls_bmd_bc_human, bc_human_results) 
cytokine_mouse_combined <- generate_combined_results(ls_bmd_cytokine_mouse, cytokine_mouse_ec5_draws_all)
cytokine_human_combined <- generate_combined_results(ls_bmd_cytokine_human, cytokine_human_ec5_draws_all)
geno_human_combined <- generate_combined_results(ls_bmd_geno_human, geno_human_results) 
ldh_human_combined <- generate_combined_results(ls_bmd_ldh_human, ldh_human_results)



# make the assay, cell type name shorter
cyto_mouse_combined <- shorten_cell_types(cyto_mouse_combined)
cyto_human_combined <- shorten_cell_types(cyto_human_combined)
bc_mouse_combined <- shorten_assay_names(bc_mouse_combined)
bc_human_combined <- shorten_assay_names(bc_human_combined)
cytokine_mouse_combined <- shorten_cell_types(cytokine_mouse_combined)
cytokine_human_combined <- shorten_cell_types(cytokine_human_combined)
geno_human_combined$assay_short = "Alkaline comet"
ldh_human_combined$assay_short = "LDH"

# ==================== 3.1 check if the data should be removed  ====================
# for mouse cyto toxicity, the largest response should larger than 10%
indices_with_results_over <- function(dflist, column = "Results", threshold = 10) {
  if (!is.list(dflist)) stop("`dflist` must be a list of data frames.")
  hits <- vapply(dflist, function(d) {
    if (!is.data.frame(d) || !(column %in% names(d))) return(FALSE)
    vals <- suppressWarnings(as.numeric(d[[column]]))
    any(vals >= threshold, na.rm = TRUE)
  }, logical(1))
  unname(which(hits))  # return unnamed integer indices
}


idx_over10_cyto_mouse <- indices_with_results_over(ls_bmd_cyto_mouse)
ls_non_trend_mouse_cyto <- setdiff(unique(cyto_mouse_combined$study), unname(idx_over10_cyto_mouse))
ls_non_trend_mouse_cyto
# 5, 17, 22, 24; 5 quite close to the EC10, 22 only 2 data points

idx_over10_cyto_human <- indices_with_results_over(ls_bmd_cyto_human)
ls_non_trend_human_cyto <- setdiff(unique(cyto_human_combined$study), unname(idx_over10_cyto_human))
ls_non_trend_human_cyto
# 1,33,42,45



# for all, check if they have over 2 order of magnitudes range
uncertain_studies <- function(df, q5="ec10_q5", q95="ec10_q95", id="study", thr=2) {
  x <- log10(df[[q5]]); y <- log10(df[[q95]])
  df[[id]][is.finite(x) & is.finite(y) & abs(x - y) >= thr]
}

ls_uncertain_mouse_cyto <- uncertain_studies(cyto_mouse_combined)
ls_uncertain_mouse_cyto
ls_uncertain_human_cyto <- uncertain_studies(cyto_human_combined)
ls_uncertain_human_cyto

ls_uncertain_mouse_bc <- uncertain_studies(bc_mouse_combined)
ls_uncertain_mouse_bc
ls_uncertain_human_bc <- uncertain_studies(bc_human_combined)
ls_uncertain_human_bc

ls_uncertain_human_geno <- uncertain_studies(geno_human_combined) 
ls_uncertain_human_geno
ls_uncertain_human_ldh  <- uncertain_studies(ldh_human_combined) 
ls_uncertain_human_ldh


ls_non_trend_mouse_cyto = c(ls_non_trend_mouse_cyto,20,21) # add also study 20, 21, in original paper it says no toxicity  https://onlinelibrary.wiley.com/doi/10.1002/smll.202000528


#======================== 4. plotting for the different EC5 value ===============
# screening out some experiment results such as for the biochemistry 
ls_non_trend_human_bc <- c(1, 4, 6, 10, 11, 12, 13, 19, 20, 21, 25, 26, 27, 29, 
                           30, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 48, 
                           49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 76)

# Add a highlight column: TRUE if the study is in ls_non_trend_human_bc
bc_human_combined <- bc_human_combined %>%
  mutate(
    highlight = ifelse(study %in% ls_non_trend_human_bc, "Non-trend", "Trend")
  )
# Filter only "Trend" studies (i.e., not in ls_non_trend_human_bc)
bc_human_combined_trend_only <- bc_human_combined %>%
  filter(!(study %in% ls_non_trend_human_bc)) %>%
  mutate(highlight = "Trend")

#remove these studies, after checking the fitted dose response curve, we found that
# the curve is almost a flat line, which is unreasonable, so we remove them
ls_further_unreasonable_human_bc = c(16,17,18,22,23,24,31,32,33,34,60,72)

bc_human_combined_trend_only <- bc_human_combined_trend_only %>%
  filter(!(study %in% ls_further_unreasonable_human_bc))

# remove the uncertain ones
bc_human_combined_trend_only <- bc_human_combined_trend_only %>%
  filter(!(study %in% ls_uncertain_human_bc))

cyto_human_combined <- cyto_human_combined %>%
  filter(!(study %in% ls_uncertain_human_cyto))%>%
  filter(!(study %in% ls_non_trend_human_cyto))

geno_human_combined <- geno_human_combined %>%
  filter(!(study %in% ls_uncertain_human_geno))
ldh_human_combined  <- ldh_human_combined %>%
  filter(!(study %in% ls_uncertain_human_ldh))


# ------------------------------- 4.1  for human --------------------------------

bc_human_combined_trend_only = convert_HD_column(bc_human_combined_trend_only)
cyto_human_combined = convert_HD_column(cyto_human_combined)
cytokine_human_combined = convert_HD_column(cytokine_human_combined)
geno_human_combined = convert_HD_column(geno_human_combined)
ldh_human_combined = convert_HD_column(ldh_human_combined)


cyto_human_combined <- cyto_human_combined %>%
  mutate(
    NP = case_when(
      str_detect(NP, regex("Graphene oxide", ignore_case = TRUE)) ~ "GO",
      TRUE ~ NP  # fallback: retain full string
    )
  )


nrow(bc_human_combined_trend_only%>%
  filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3")))
nrow(cyto_human_combined%>%
       filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3")))
nrow(cytokine_human_combined%>%
       filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3")))
nrow(geno_human_combined%>%
       filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3")))
nrow(ldh_human_combined%>%
       filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3")))


cyto_human_combined <- cyto_human_combined %>%
  mutate(general_experiment_type ="Cytotoxicity",
          experiment_type = "Cytotoxicity",
         species = "human")

bc_human_combined_trend_only <- bc_human_combined_trend_only %>%
  mutate(general_experiment_type ="Biochemical",
          experiment_type = "Biochemical",
         species = "human")

cytokine_human_combined <- cytokine_human_combined %>%
  mutate(general_experiment_type ="Cytokine",
         experiment_type = "Cytokine",
         species = "human")

cytokine_human_combined <- cytokine_human_combined %>%
  mutate(
    ec10_median = as.numeric(na_if(min_signif_conc, "-"))
  )

geno_human_combined <- geno_human_combined %>%
  mutate(general_experiment_type ="Genotoxicity",
        experiment_type = "Genotoxicity",
         species = "human")

ldh_human_combined <- ldh_human_combined %>%
  mutate(general_experiment_type ="Cytotoxicity",
      experiment_type = "LDH",
         species = "human")

combined_human_df <- bind_rows(cyto_human_combined, bc_human_combined_trend_only, 
                               cytokine_human_combined,geno_human_combined,
                               ldh_human_combined)

combined_human_df <- combined_human_df %>%
  mutate(
    NP = case_when(
      #str_detect(NP, regex("Graphene quantum dots", ignore_case = TRUE)) ~ "GQD",
      str_detect(NP, regex("Graphene oxide", ignore_case = TRUE)) ~ "GO",
      TRUE ~ NP  # fallback: retain full string
    )
  )
save(combined_human_df, 
     file = file.path("~/work/code/codo_v2/R/results_undissolved/human/stan_results/combined_human_ec_results.RData"))
#load(file = file.path("~/work/code/codo_v2/R/results_undissolved/human/stan_results/combined_human_ec_results.RData"))


# remove the ratio when the range for EC10 is too large if it covers 3 magnitude
combined_human_df_filtered <- combined_human_df %>%
  filter(
    experiment_type == "Cytokine" |
      (experiment_type != "Cytokine" & abs(log10(ec10_q5) - log10(ec10_median)) <= 3)
  )
combined_human_df_filtered
save(combined_human_df_filtered, 
     file = file.path("~/work/code/codo_v2/R/results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData"))


# ================================== 4.2  for mouse ================================== 
bc_mouse_combined = convert_HD_column(bc_mouse_combined)
cyto_mouse_combined = convert_HD_column(cyto_mouse_combined)
cytokine_mouse_combined = convert_HD_column(cytokine_mouse_combined)


cyto_mouse_combined = cyto_mouse_combined %>%
  filter(!(study %in% ls_uncertain_mouse_cyto))%>%
  filter(!(study %in% ls_non_trend_mouse_cyto))


bc_mouse_combined = bc_mouse_combined %>%
  filter(!(study %in% ls_uncertain_mouse_bc))



cyto_mouse_combined <- cyto_mouse_combined %>%
  mutate(general_experiment_type = "Cytotoxicity",
    experiment_type = "Cytotoxicity",
         species = "mouse")

bc_mouse_combined <- bc_mouse_combined %>%
  mutate(general_experiment_type = "Biochemical",
    experiment_type = "Biochemical",
         species = "mouse")

cytokine_mouse_combined <- cytokine_mouse_combined %>%
  mutate(general_experiment_type = "Cytokine",
    experiment_type = "Cytokine",
         species = "mouse")

cytokine_mouse_combined <- cytokine_mouse_combined %>%
  mutate(
    ec10_median = as.numeric(na_if(min_signif_conc, "-"))
  )%>%filter(!is.na(ec10_median))


combined_mouse_df <- bind_rows(cyto_mouse_combined, bc_mouse_combined,cytokine_mouse_combined)


combined_mouse_df <- combined_mouse_df %>%
  mutate(
    NP = case_when(
      str_detect(NP, regex("Amorphous SiO2", ignore_case = TRUE)) ~ "SiO2",
      str_detect(NP, regex("Graphene oxide", ignore_case = TRUE)) ~ "GO",
      TRUE ~ NP  # fallback: retain full string
    )
  )
#save(combined_mouse_df, 
#     file = file.path("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_results.RData"))

combined_mouse_df_filtered <- combined_mouse_df %>%
  filter(
    experiment_type == "Cytokine" |
      (experiment_type != "Cytokine" & abs(log10(ec10_median) - log10(ec10_q5)) <= 3)
  )
save(combined_mouse_df_filtered, 
     file = file.path("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_filtered_results.RData"))

