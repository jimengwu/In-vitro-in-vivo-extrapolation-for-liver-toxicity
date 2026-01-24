# this function will first read the collected in vitro dose response data, and then clean
# the data set, making them the right format for next processing step 
# (single one by one experiment record)

# Load necessary libraries
library(readxl)  # for reading Excel files
library(dplyr)   # for data manipulation
library(proast71.1)


#---------------------A. calculate the BMD for the human/mouse liver in vitro data -----------


# -------------- 0. Read the Excel file ------------
# human liver cell data
liver_results <- read_excel("/Users/wuji/work/code/codo_v2/calculation_in_vitro/human/processed_human_liver_0716.xlsx")
outputfile = "results_undissolved/human/ls_sub_case_human.RData"
proast_folder = "/Users/mmm/work/code/codo_v2/R/results_undissolved/human/bmd_results_undissolved/"

# mouse liver cell data 
liver_results <- read_excel("/Users/wuji/work/code/codo_v2/calculation_in_vitro/mouse/processed_mouse.xlsx")
outputfile = "results_undissolved/mouse/ls_sub_case_mouse.RData"
proast_folder = "/Users/mmm/work/code/codo_v2/R/results_undissolved/mouse/bmd_results_undissolved/"


#------------- 0.1 unify the concentration unit --------------

#------------unify the concentration unit into ug/cm3--------
# Assuming liver_results is your data frame in R

# Remove rows with NA values in 'In vitro concentration' or 'In vitro concentration unit'
liver_results <- liver_results[!is.na(liver_results$`In vitro concentration unit`), ]

# Delete rows where "In vitro concentration unit" is "nM"
liver_results <- liver_results[liver_results$`In vitro concentration unit` != "nM", ]

# Check unique values in "In vitro concentration unit"
unique(liver_results$`In vitro concentration unit`)

# Convert "mg/ml" to "μg/cm3" (1 mg/ml = 1000 μg/cm3)
liver_results$`In vitro concentration` <- as.numeric(liver_results$`In vitro concentration`)
liver_results$`In vitro concentration`[liver_results$`In vitro concentration unit` == "mg/ml"] <- 
  liver_results$`In vitro concentration`[liver_results$`In vitro concentration unit` == "mg/ml"] * 10^3

liver_results$`In vitro concentration`[liver_results$`In vitro concentration unit` == "mg/cm3"] <- 
  liver_results$`In vitro concentration`[liver_results$`In vitro concentration unit` == "mg/cm3"] * 10^3



# For "μg/cm2" and "96-well culture plate" (area = 0.32 cm2, volume = 0.1 ml)
liver_results$`In vitro concentration`[
  liver_results$`In vitro concentration unit` == "μg/cm2" & liver_results$`Type of plate` == "96-well culture plate"
] <- liver_results$`In vitro concentration`[
  liver_results$`In vitro concentration unit` == "μg/cm2" & liver_results$`Type of plate` == "96-well culture plate"
] * (0.32 / 0.1)

# For "μg/cm2" and "6-well culture plate" (area = 9.6 cm2, volume = 1 ml)
liver_results$`In vitro concentration`[
  liver_results$`In vitro concentration unit` == "μg/cm2" & liver_results$`Type of plate` == "6-well culture plate"
] <- liver_results$`In vitro concentration`[
  liver_results$`In vitro concentration unit` == "μg/cm2" & liver_results$`Type of plate` == "6-well culture plate"
] * (9.6 / 1)

# Change all "In vitro concentration unit" to "μg/cm3"
#liver_results$`In vitro concentration unit` <- "μg/cm3"

#------------ 0.2. get only numeric results -------



liver_results <- liver_results[liver_results$`Raw Results` != "-", ]

# Remove the " ±" and everything after it, then convert to numeric
liver_results$`Raw Results` <- as.numeric(gsub(" ±.*", "", liver_results$`Raw Results`))
liver_results$`In vitro concentration` <- as.numeric(liver_results$`In vitro concentration`)

# Check if 'Tested assay' contains the string "viability" and update 'Raw Results'
liver_results$Results <- liver_results$`Raw Results`

# Step 1: get "viability" rows
viab_idx <- grepl("viability", liver_results$`Tested assay`, ignore.case = TRUE) & grepl("%", liver_results$`Tested assay`, ignore.case = TRUE)

liver_results$Results[viab_idx] <- 100 - liver_results$Results[viab_idx]

# Step 2: Remove rows where Results <= 0 AND Tested assay contains "viability" 
# what was deleted here are all the cytotoxicity results, which means there are no clear trend in the toxicity 
# the negstive means more cells are aviable given the np concentration
# for biochemistry, minimal value is needed to be normalized the response

liver_results$Results[
  liver_results$Results == 0 & grepl("Biochemical", liver_results$`Tested assay`, ignore.case = TRUE)
] <- 1e-6
liver_results <- liver_results[!(liver_results$Results <= 0), ]  # for cytotoxicty no zero reponse 

# Remove rows where both 'Raw Results' and 'In vitro concentration' are zero since the proast will have problems on this

liver_results <- liver_results[!(liver_results$`Results` == 0 & liver_results$`In vitro concentration` == 0), ]

# Check the updated liver_results to ensure the rows are removed
print(nrow(liver_results))  # Print the number of remaining rows





# -------------- 1. split the raw database into the format as list with single experiment record one by one ------------
# Initialize the list to hold subsets
ls_sub_case <- list()  
ls_zero_conc_case <- list()
library(dplyr)

# for mouse data
liver_results <- liver_results %>%
  mutate(`Specific surface area (m2/g)` = ifelse(is.na(`Specific surface area (m2/g)`), "no_value", `Specific surface area (m2/g)`))

unique_cores <- liver_results %>%
      distinct(`Substance name`) %>%
       filter(!is.na(`Substance name`)) %>%
       arrange(`Substance name`) %>%
       pull()

allowed_cores <- c(
  "Au",
  "Graphene oxide",
  "Graphene quantum dots (GQD)",
  "SiO2",
  "TiO2",
  "rGO",
  "Amorphous SiO2",
  "Fe2O3"
)
# Loop over unique Study_IDs to make the single data frame into a list of each experiment
for (i in unique(liver_results$Study_ID)) {
  

   # Filter by Study_ID
  case <- liver_results[liver_results$Study_ID == i, ]

  # remove rows where all values are nan
  case <- case[!apply(case, 1, function(row) all(is.na(row) | is.nan(row))), ]
  
  
  split_cases <- split(case, list(case$Study_ID,case$`Substance core`, 
                                  case$`Tested assay`,case$`Cell type`,
                                  case$`Average size (nm)`,
                                 case$`Specific surface area (m2/g)`, # this is for mouse in vitro data, in study 7 two Tio2 has the same diameter but different ssa
                                  case$`In vitro media`,case$`In vitro exposure time[h]`))
  
  # Convert each element in split_cases to a standard data frame
  split_cases <- lapply(split_cases, as.data.frame)
  
  split_cases <- Filter(function(df) nrow(df) > 0, lapply(split_cases, as.data.frame)) # remove empty dataframe
  
  # Assuming split_cases is your list of data frames
  single_zero_conc <- lapply(split_cases, function(df) {
    if (nrow(df) == 1 && df$`In vitro concentration` == 0) return(df)
    else return(NULL)
  })
  
  # Remove NULLs (i.e., data frames that didn't match the condition)
  single_zero_conc <- Filter(Negate(is.null), single_zero_conc)
  
  ls_zero_conc_case <- append(single_zero_conc, ls_zero_conc_case)
  
  # Filter out data frames that are NOT both one-row and have zero concentration
  processed_sub_cases <- Filter(function(df) {
    !(nrow(df) == 1 && df$`In vitro concentration` == 0)
  }, split_cases)
  
  # Append the non-null cases to ls_sub_case
  if (length(processed_sub_cases) > 0) {
    ls_sub_case <- append(ls_sub_case, processed_sub_cases)
    
  rm(processed_sub_cases)
  }

}



# create the file name for each experiment for the subsequent dose response analysis results saving 
for (i in 1:length(ls_sub_case)){
 
  temp = ls_sub_case[[i]]

  combined_unique_values <- c(unique(temp$Study_ID),unique(temp$`Substance core`),
                                     unique(temp$`Tested assay`),unique(temp$`Cell type`),
                                     unique(temp$`Average size (nm)`),
                                     #unique(temp$`Specific surface area (m2/g)`), # this is for mouse in vitro data, in study 7 two Tio2 has the same diameter but different ssa
                                     unique(temp$`In vitro media`),
                                     unique(temp$`In vitro exposure time[h]`))
  
  combined_string <- paste(combined_unique_values, collapse = "_")
  
  safe_combined_string <- gsub("[^[:alnum:]_-]", "_", combined_string)
  #print(combined_string)
  combined_file = paste0(proast_folder,safe_combined_string,".RData")
  ls_sub_case[[i]]$file = combined_file
  print(combined_file)

}


keep <- vapply(ls_sub_case, function(temp) {
  core <- unique(na.omit(temp$`Substance name`))
  length(core) == 1 && core %in% allowed_cores
}, logical(1))

ls_sub_case <- ls_sub_case[keep]
# create parent directory if missing
dir.create(dirname(outputfile), recursive = TRUE, showWarnings = FALSE)

save(ls_sub_case, file = outputfile)

