# ==============================================================================
# TITLE:   Collection and visualisation of the BMD10 results across all
#          in vitro experiments (mouse and human)
#
# PURPOSE: Reads the per-experiment EC10 (BMD10) output files produced by the
#          endpoint-specific fitting scripts (cytotoxicity, metabolic cell
#          stress, genotoxicity) together with the LOEC-based cytokine
#          thresholds, and merges them into a single harmonised results table.
#          Each record is annotated with its experiment identity (nanoparticle
#          core, cell line, assay, assay category) and species, and the median
#          estimate is reported with its credible interval.
#          The script then generates the summary figures: distribution of EC10
#          values (ug/mL) across nanoparticle types, plotted separately for
#          human and mouse cell systems, as box-whisker plots with individual
#          experiments shown as points and coloured by assay category.
#
# INPUT:   <EC10 / BMD10 result files from scripts 02-04, one per endpoint class>
# OUTPUT:  <combined EC10 summary table (mouse + human)>
#          <human EC10 distribution figure; mouse EC10 distribution figure>
#          <per-experiment credible interval figures, faceted by nanoparticle
#           and cell type>
#
# NOTE:    Cytokine values are LOEC-based surrogates rather than fitted EC10
#          estimates; they are shown for completeness and flagged as such.
#
# AUTHOR:  Jimeng Wu          CREATED: 2026-04          
# ==============================================================================


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
#ls_non_trend_mouse_cyto <- setdiff(ls_non_trend_mouse_cyto, c(5L))

idx_over10_cyto_human <- indices_with_results_over(ls_bmd_cyto_human)
ls_non_trend_human_cyto <- setdiff(unique(cyto_human_combined$study), unname(idx_over10_cyto_human))
ls_non_trend_human_cyto
#ls_non_trend_human_cyto <- setdiff(ls_non_trend_human_cyto, c(33L, 45L))
# 1,33,42,45

# but if the generated dose response curve has a low uncertainty range and the highest response close to 10, 
# we would still keep that dataset and corresponding ecx results
# for human: study 53, 76
# for mouse: study 5, 19, 38; 23, 36 only two data points.

#ls_non_trend_human_cyto <- setdiff(ls_non_trend_human_cyto, c(53L, 76L))
#ls_non_trend_mouse_cyto <- setdiff(ls_non_trend_mouse_cyto, c(5L, 19L, 38L))



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

#ls_uncertain_human_cyto<- setdiff(ls_uncertain_human_cyto, c(40L, 41L))

ls_non_trend_mouse_cyto = c(ls_non_trend_mouse_cyto,20,21) # add also study 20, 21, in original paper it says no toxicity  https://onlinelibrary.wiley.com/doi/10.1002/smll.202000528

# undissolve
# cyto mouse uncertainty: 16 23
# cyto mouse smaller than 10: 5 17 22 24
# cyto human: 19 23 24 25 26 29 30 31 32 40 41 42 # 40,41 close to EC10
# BC HUMAN: 34 72 74

# all
#ls_uncertain_human_cyto<- setdiff(ls_uncertain_human_cyto, c(68L, 69L))
# uncertainty range has the problem of responses accumulated together.
# cyto human: study 38, 42, 43, 44, 45, 48, 49, 50, 51, 71; 
# study 68,69 close to the measured concentrations, maybe could be kept

# geno human, study: 14,15 seems has the same problem as in the cytotoxicity, which they only
# have concentrations that showing larger response compared to benchmark
# bc human, study 203, 193, 192, 157, 151 also have the same problem in the cytotoxicity, which they have 
# only concentrations that showing larger response compared to benchmark



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
ls_further_unreasonable_human_bc = c(16,17,18,22,23,24,31,32,33,34,60,72,47)

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
      #str_detect(NP, regex("Graphene quantum dots", ignore_case = TRUE)) ~ "GQD",
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
  )%>%filter(!is.na(ec10_median))


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

# ============================ 4.1.B read the processed human dose response curve data ==============================
# load the filtered object combined_human_df_filtered
load(file = file.path("~/work/code/codo_v2/R/results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData"))

# check the filtered id 

# Studies in the full set but not in the filtered set

all_pairs      <- unique(paste(combined_human_df$study, combined_human_df$assay, sep = "_"))
filtered_pairs <- unique(paste(combined_human_df_filtered$study, combined_human_df_filtered$assay, sep = "_"))

missing_pairs <- setdiff(all_pairs, filtered_pairs)
# Optionally, view them
missing_pairs


library(ggplot2)
library(forcats)
#plot_ec10_human = plot_ecx_human(combined_human_df_filtered,"ec10")  # For EC10
#plot_ec5_human = plot_ecx_human(combined_human_df_filtered,"ec5")  # For EC10
folder = "results_undissolved/plots/before_extrapolation/"
#ggsave(paste0(folder,"ec5_human", ".pdf"), plot = plot_ec5_human, width = 8, height = 9)
#ggsave(paste0(folder,"ec10_human", ".pdf"), plot = plot_ec10_human, width = 8, height = 9)


folder = "results_undissolved/plots/before_extrapolation/"

combined_human_df_filtered_undissolved = combined_human_df_filtered %>%
  filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3"))

subset_df <- combined_human_df_filtered_undissolved[
  combined_human_df_filtered_undissolved$NP == "SiO2" &
    combined_human_df_filtered_undissolved$Cell_type == "L-02" &
    combined_human_df_filtered_undissolved$general_experiment_type == "Cytotoxicity",
]

combined_human_df_filtered_undissolved$general_experiment_type[
  combined_human_df_filtered_undissolved$general_experiment_type == "Biochemical"
] <- "Metabolic Cell Stress"

combined_human_df_filtered_undissolved$experiment_type[
  combined_human_df_filtered_undissolved$experiment_type == "Biochemical"
] <- "Metabolic Cell Stress"

plot_ec10_human_undisolved = plot_ecx_human(combined_human_df_filtered_undissolved,"ec10")  # For EC10

#ggsave(paste0(folder,"ec5_human_undisolved", ".pdf"), plot = plot_ec5_human_undisolved, width = 8, height = 9)
#ggsave(paste0(folder,"ec10_human_undisolved", ".pdf"), 
#       plot = plot_ec10_human_undisolved, width = 7, height = 5,device = cairo_pdf)

plot_ec10_human_undisolved_main = plot_ecx_human_main(combined_human_df_filtered_undissolved,"ec10")  # For EC10
plot_ec10_human_undisolved_main
ggsave(paste0(folder,"ec10_human_undisolved_combined", ".pdf"), 
       plot = plot_ec10_human_undisolved_main, width = 6, height = 4,device = cairo_pdf)

plot_ec10_human_undisolved_si = plot_ecx_human_si(combined_human_df_filtered_undissolved,"ec10")  # For EC10
plot_ec10_human_undisolved_si
ggsave(paste0(folder,"ec10_human_undisolved_separated", ".pdf"), 
       plot = plot_ec10_human_undisolved_si, width = 9, height = 16,device = cairo_pdf)


df_subset <- combined_human_df_filtered_undissolved %>%
  filter(NP == "SiO2", grepl("cytotoxicity", assay, ignore.case = TRUE))

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
# ======================================== 4.2.B  read data for mouse ========================================
# load the filtered object combined_mouse_df_filtered
load(file = file.path("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_results.RData"))
# check the filtered id ----------------------------


all_pairs_mouse      <- unique(paste(combined_mouse_df$study, combined_mouse_df$assay, sep = "_"))
filtered_pairs_mouse <- unique(paste(combined_mouse_df_filtered$study, combined_mouse_df_filtered$assay, sep = "_"))

missing_pairs_mouse <- setdiff(all_pairs_mouse, filtered_pairs_mouse)
# Optionally, view them
missing_pairs_mouse


# check the filtered id ---------


#plot_ec10_mouse = plot_ecx_human(combined_mouse_df_filtered,"ec10_median")  # For EC10
#plot_ec5_mouse = plot_ecx_human(combined_mouse_df_filtered,"ec5_median")  # For EC5
#folder = "results/plots/before_extrapolation/"
#ggsave(paste0(folder,"ec5_mouse", ".pdf"), plot = plot_ec5_mouse, width = 8, height = 9)
#ggsave(paste0(folder,"ec10_mouse", ".pdf"), plot = plot_ec10_mouse, width = 8, height = 9)


combined_mouse_df_filtered_undissolved = combined_mouse_df %>%
  filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3"))

combined_mouse_df_filtered_undissolved$general_experiment_type[
  combined_mouse_df_filtered_undissolved$general_experiment_type == "Biochemical"
] <- "Metabolic Cell Stress"

combined_mouse_df_filtered_undissolved$experiment_type[
  combined_mouse_df_filtered_undissolved$experiment_type == "Biochemical"
] <- "Metabolic Cell Stress"

plot_ec10_mouse_undisolved = plot_ecx_human(combined_mouse_df_filtered_undissolved,"ec10")  # For EC10
folder = "results_undissolved/plots/before_extrapolation/"
#ggsave(paste0(folder,"ec10_mouse_undisolved", ".pdf"), plot = plot_ec10_mouse_undisolved, 
#       width = 7, height = 5,device = cairo_pdf)

plot_ec10_mouse_undisolved_main = plot_ecx_human_main(combined_mouse_df_filtered_undissolved,"ec10")  # For EC10
plot_ec10_mouse_undisolved_main

ggsave(paste0(folder,"ec10_mouse_undisolved_combined", ".pdf"), 
       plot = plot_ec10_mouse_undisolved_main, width = 6, height = 4,device = cairo_pdf)

plot_ec10_mouse_undisolved_si = plot_ecx_human_si(combined_mouse_df_filtered_undissolved,"ec10")  # For EC10
plot_ec10_mouse_undisolved_si
ggsave(paste0(folder,"ec10_mouse_undisolved_separated", ".pdf"), 
       plot = plot_ec10_mouse_undisolved_si, width = 10, height = 8,device = cairo_pdf)




# ===================== 5. analyze the SiO2 size relationship =====================
# load the filtered object combined_human_df_filtered
load(file = file.path("~/work/code/codo_v2/R/results_undissolved/human/stan_results/combined_human_ec_filtered_results.RData"))
load(file = file.path("~/work/code/codo_v2/R/results_undissolved/mouse/stan_results/combined_mouse_ec_filtered_results.RData"))
# filter the undissolved NPs
combined_human_df_filtered_undissolved = combined_human_df_filtered %>%
  filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3"))
combined_mouse_df_filtered_undissolved = combined_mouse_df_filtered %>%
  filter(NP %in% c("Au", "TiO2", "SiO2","GO","rGO","Fe2O3"))

combo_counts <- combined_human_df_filtered_undissolved %>%
  count(NP, Cell_type,general_experiment_type, name = "n") %>%
  arrange(desc(n))

df_freq = combined_human_df_filtered_undissolved %>%
  count(NP, general_experiment_type, Cell_type, name = "n")


#filter the TiO2
combined_human_df_NP = combined_human_df_filtered_undissolved %>%
  filter(NP == "TiO2") %>%
  filter("Cytotoxicity" == general_experiment_type) %>%
  filter(Cell_type == "HepG2") %>%
  filter(!is.na(Average_size_nm)) %>%
  mutate(Average_size_nm = as.numeric(Average_size_nm))

combined_human_df_NP = convert_HD_column(combined_human_df_NP)

df_filtered <- combined_human_df_NP %>%
  select(study, Paper_ID,NP_core, exposure_time, experiment_type, assay,
         Average_size_nm, HD,ec10_median)

library(dplyr)
library(ggplot2)
library(patchwork)

# ---- regression model ----
model <- lm(log10(ec10_median) ~ Average_size_nm, data = combined_human_df_NP)
r2 <- summary(model)$r.squared
pval <- summary(model)$coefficients[2,4]

# ---- Scatter with regression ----
p1 <- ggplot(combined_human_df_NP,
             aes(x = HD, y = ec10_median)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, fill = "#D5D1D1") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  scale_y_log10(labels = scales::label_scientific(),
                limits = c(1e-1, 1e3)) +
  labs(
    x = "Hydrodynamic diameter (nm)",
    y = expression(bold("Median"~ bolditalic("in vitro") ~BMD[10]* " (µg/cm"^3*")"))
  )+
  #annotate("text",
  #         x = max(combined_human_df_Sio2$Average_size_nm) * 0.7,
  #         y = max(combined_human_df_Sio2$ec10_median, na.rm = TRUE),
  #         label = paste0("R² = ", round(r2, 2), 
  #                        "\n", "p = ", signif(pval, 2)),
  #         hjust = 0, vjust = 1, size = 5) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16,face = "bold",color = "black"),   # axis titles
    axis.text  = element_text(size = 14,color = "black"),                  # axis tick labels
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)  # <- panel border
  )

# --- Boxplot with adjusted bins and borders ---
p2 <- combined_human_df_NP %>%
  mutate(size_bin = cut(HD,
                        breaks = c(0, 100, 200,1000),
                        labels = c("≤100 nm", "101–200 nm", "201–1000 nm"))) %>%
  ggplot(aes(x = size_bin, y = ec10_median, fill = size_bin)) +
  # violin
  geom_violin(alpha = 0.6, trim = FALSE, color = "black", size = 0.5) +
  # boxplot inside
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
               color = "black", size = 0.4) +
  # jittered points
  geom_jitter(width = 0.15, alpha = 0.7, size = 2, shape = 21) +
  scale_fill_manual(values = c("≤100 nm"   = "#BFDFD2",
                               "101–200 nm" = "#ECB66C",
                               "201–1000 nm"= "#ED8D5A")) +
  
  #scale_fill_manual(values = c("≤30 nm"   = "#B6B3D6",
  #                             "31–100 nm" = "#F8B2A2",
  #                             "101–200 nm"= "#E9587A")) +
  scale_y_log10(labels = scales::label_scientific(),
                limits = c(1e-1, 1e3)) +
  labs(
    x = "Hydrodynamic diameter range (nm)",
    y = expression(bold("Median"~ bolditalic("in vitro") ~ BMD[10]* " (µg/cm"^3*")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16,face = "bold",color = "black"),   # axis titles
    axis.text  = element_text(size = 14,color = "black"),                  # axis tick labels
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)  # <- panel border
  )

# ---- Combine plots side by side ----
final_plot <- p2 + p1 +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"))

final_plot

ggsave("results_undissolved/plots/before_extrapolation/TiO2_size_vs_EC10_HepG2.pdf",
       plot = final_plot, width = 11, height = 5,
       device = cairo_pdf)

color_list = c("#BFDFD2","#51999f","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A")

