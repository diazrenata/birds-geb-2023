if(FALSE) {
  
  install.packages("dplyr")
  install.packages("tidyr")
  install.packages("vegan")
  install.packages("here")
  install.packages("mclust")
  install.packages("MuMIn")
  install.packages("R.utils")
}

library(dplyr)

R.utils::sourceDirectory(here::here("R"))

datasets_files <- list.files(here::here("analysis", "bbs_data"), pattern = ".RDS")

# select datasets to work with
working_datasets <- read.csv(here::here("analysis", "supporting_data","eightypercent_coverage_1988_2018.csv"))

working_datasets$filename <- paste0(working_datasets$matssname, ".RDS")

datasets_to_use <- datasets_files[ datasets_files %in% working_datasets$filename]

datasets_to_use <- datasets_to_use[1:200]

# lapply to load datasets

datasets_loaded <- lapply(datasets_to_use, FUN = function(x) readRDS(here::here("analysis", "bbs_data", x)))

set.seed(1989)

# Generate and calculate estimated annual total biomass, energy use, and abundance per route
annual_svs <- lapply(datasets_loaded, FUN = get_annual_state_variables)

# Combine to a single data frame and save
all_sims <- dplyr::bind_rows(annual_svs)
saveRDS(all_sims, file = here::here("analysis", "results", "all_sims.RDS"))

# Fit GLMS, calculate AIC values, and calculate model predictions for biomass
# Save AICs and model predictions

glms_b <- lapply(annual_svs, FUN = fit_trend_models_biomass)

aics_b <- lapply(glms_b, FUN = eval_trend_models)
all_aics_b <- dplyr::bind_rows(aics_b)
saveRDS(all_aics_b, file = here::here("analysis", "results", "all_aics_b.RDS"))

preds_b <- lapply(glms_b, FUN = all_models_predicted_change)
all_preds_b <- dplyr::bind_rows(preds_b)
saveRDS(all_preds_b, file = here::here("analysis", "results", "all_preds_b.RDS"))

# Do the same for energy use

glms_e <- lapply(annual_svs, FUN = fit_trend_models_energy)

aics_e <- lapply(glms_e, FUN = eval_trend_models)
all_aics_e <- dplyr::bind_rows(aics_e)
saveRDS(all_aics_e, file = here::here("analysis", "results", "all_aics_e.RDS"))

preds_e <- lapply(glms_e, FUN = all_models_predicted_change)
all_preds_e <- dplyr::bind_rows(preds_e)
saveRDS(all_preds_e, file = here::here("analysis", "results", "all_preds_e.RDS"))

# Compare first and last five-year periods for each dataset

cs_compares <- lapply(datasets_loaded, FUN = compare_community_structure)
all_cs_compares <- dplyr::bind_rows(cs_compares)
saveRDS(all_cs_compares, file = here::here("analysis", "results", "all_cs_compares.RDS"))
