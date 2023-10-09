#' Slice data
#'
#' @param dataset a MATSS-formatted dataset
#' @param keep_years which years to keep
#'
#' @return dataset subsetted to the years in keep_years
#' @export
#'
pull_focal_years <- function(dataset, keep_years = c(1989:2018)) {

  keep_rows <- which(dataset$covariates$year %in% keep_years)

  dat_out <- dataset

  dat_out$abundance <- dat_out$abundance[keep_rows, ]
  dat_out$covariates <- dat_out$covariates[keep_rows,]

  dat_out
}

#' Null model of no size change
#'
#' Resamples species' abundances according to their mean relative abundance across the entire timeseries.
#'
#' This is the null model to characterize "no change in the size structure" by keeping species composition fixed to the long-term average across the timeseries.
#'
#' @param dataset a MATSS dataset
#' @param resp_seed a seed
#'
#' @return df
#' @export
relabund_null_model <- function(dataset, resp_seed = 1989) {

  # Compute each species' relative abundance (proportion of total abundance in that year) in each year
  spRelAbunds <- dataset$abundance/ rowSums(dataset$abundance)

  # Compute each species' average relative abundance across all years
  spRelAbunds <- colSums(spRelAbunds) / nrow(spRelAbunds)

  # Function to draw N individuals (N being the number of individuals seen in a given year) with probabilities weighted according to species' mean relative abundance across all eyears
  draw_year <- function(year_row, spRelAbunds) {
    as.data.frame(t(rmultinom(1, size  = sum(year_row), prob = spRelAbunds)))
  }

  set.seed(resp_seed)

  # Draw species' abundances
  newAbund <- apply(dataset$abundance, MARGIN = 1, FUN = draw_year, spRelAbunds = spRelAbunds)
  set.seed(NULL)

  dataset$abundance <- dplyr::bind_rows(newAbund)

  dataset
}

#' Calculate state variables for each year off an ISD
#'
#' @param isd list with $isd, $metadata. $metadata is MATSS metadata. $isd is a dataframe with columns id (spID), mass, year, isd_seed
#'
#' @return dataframe with state variables off ISD plus location data
#' @export
#'
#' @importFrom dplyr mutate group_by summarize n ungroup rename bind_cols
calc_sv <- function(isd) {

  sv_isd <- isd$isd %>%
    dplyr::mutate(energy = estimate_b(mass)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(abundance = dplyr::n(),
                     total_biomass = sum(mass),
                     total_energy = sum(energy)) %>%
    dplyr::ungroup()%>%
    dplyr::rename(timeperiod = year)

  dplyr::bind_cols(sv_isd, isd$metadata$location)
}

#' Get annual state variables
#'
#' High level wrapper to get annual state variables for observed and null model (abundance-driven) dynamics.
#'
#' @param dat MATSS dataset
#'
#' @return df
#' @export
#' @importFrom dplyr mutate bind_rows
get_annual_state_variables <- function(dat) {

  dat_sliced <- pull_focal_years(dat) # pull years we want

  dat_resp <- relabund_null_model(dat_sliced, resp_seed = 1989) # null model

  dat_isd <- simulate_isd_ts(dat_sliced, isd_seed = 1989) # isd for observed

  dat_resp_isd <- simulate_isd_ts(dat_resp, isd_seed = 1989) # isd for null model

  # annual state variables for observed
  real_sv <- calc_sv(dat_isd) %>%
    dplyr::mutate(source = "real")

  # annual state variables for null model
  sim_sv <- calc_sv(dat_resp_isd) %>%
    dplyr::mutate(source = "sim")

  # combine annual state variables
  both_sv <- dplyr::bind_rows(real_sv, sim_sv)  %>%
    dplyr::mutate(fsource = as.factor(source),
           matssname = paste0("bbs_rtrg_", dat$metadata$location$route, "_", dat$metadata$location$statenum),
           simtype = "actual")

  both_sv
}

#' Fit models - biomass
#'
#' Fit models corresponding to syndromes of "no trend", "trend", "decoupled trend". Fits using both Gaussian and Gamma (log link) generalized linear models.
#'
#' @param sims result of get_annual_state_variables
#'
#' @return df
#' @export
#'
fit_trend_models_biomass <- function(sims) {


  models <- list(

    gaussian_glm_full = glm(total_biomass ~ timeperiod * source, data= sims),
    gaussian_glm_noint = glm(total_biomass ~ timeperiod + source, data= sims),
    gaussian_glm_nos = glm(total_biomass ~ timeperiod, data= sims),
    gaussian_glm_notime = glm(total_biomass ~ 1, data = sims),


    gamma_glm_full = glm(total_biomass ~ timeperiod * source, family = Gamma(link = "log"), data= sims),
    gamma_glm_noint = glm(total_biomass ~ timeperiod + source, family = Gamma(link = "log"), data= sims),
    gamma_glm_nos = glm(total_biomass ~ timeperiod, family = Gamma(link = "log"), data= sims),
    gamma_glm_notime = glm(total_biomass ~ 1, family = Gamma(link = "log"), data = sims)

  )

  models

}

#' Fit models - energy
#'
#' Fit models corresponding to syndromes of "no trend", "trend", "decoupled trend". (also fits a separate intercept model but in practice it never gets used). Fits using both Gaussian and Gamma (log link) generalized linear models.
#'
#' @param sims result of get_annual_state_variables
#'
#' @return list of 8 models
#' @export
#'
fit_trend_models_energy <- function(sims) {


  models <- list(

    gaussian_glm_full = glm(total_energy ~ timeperiod * source, data= sims),
    gaussian_glm_noint = glm(total_energy ~ timeperiod + source, data= sims),
    gaussian_glm_nos = glm(total_energy ~ timeperiod, data= sims),
    gaussian_glm_notime = glm(total_energy ~ 1, data = sims),


    gamma_glm_full = glm(total_energy ~ timeperiod * source, family = Gamma(link = "log"), data= sims),
    gamma_glm_noint = glm(total_energy ~ timeperiod + source, family = Gamma(link = "log"), data= sims),
    gamma_glm_nos = glm(total_energy ~ timeperiod, family = Gamma(link = "log"), data= sims),
    gamma_glm_notime = glm(total_energy ~ 1, family = Gamma(link = "log"), data = sims)

  )

  models

}


#' Model comparison metrics
#'
#' Compute p values and AIC values comparing models fit with different terms. P values could be used in a kind of stepwise selection context, but I prefer to use AIC because you can compare all 4 together instead of the sequential P values.
#'
#' @param some_models result of fit_trend_models_energy or fit_trend_models_biomass
#'
#' @return dataframe of model name, aic, p values for sequential, model complexity (1-4 with 4 high)
#' @export
#'
#' @importFrom dplyr bind_rows mutate
eval_trend_models <- function(some_models) {

  p_interaction <- anova(some_models$gaussian_glm_full, some_models$gaussian_glm_noint, test = "F")[2,6] # p value comparing year * source to year + source
  p_sourceintercept <- anova(some_models$gaussian_glm_noint, some_models$gaussian_glm_nos, test = "F")[2,6] # p value comparing year + source to year
  p_slope <- anova(some_models$gaussian_glm_nos,some_models$gaussian_glm_notime, test = "F")[2,6]  # p value comparing year to intercept-only

  gaussian_aics <- dplyr::bind_rows(lapply(some_models[1:4], pull_model_aic)) %>% # AIC for each model
    dplyr::mutate(model_p = c(p_interaction, p_sourceintercept, p_slope, NA),
           modelcomplexity = c(4,3,2,1)) # assign a model complexity for sorting purposes


  p_gamma_interaction <- anova(some_models$gamma_glm_full, some_models$gamma_glm_noint, test = "F")[2,6]# p value comparing year * source to year + source
  p_gamma_sourceintercept <- anova(some_models$gamma_glm_noint, some_models$gamma_glm_nos, test = "F")[2,6]# p value comparing year + source to year
  p_gamma_slope <- anova(some_models$gamma_glm_nos,some_models$gamma_glm_notime, test = "F")[2,6]# p value comparing year to intercept-only


  gamma_aics <- bind_rows(lapply(some_models[5:8], pull_model_aic)) %>% # AIC for each model
    mutate(model_p = c(p_gamma_interaction, p_gamma_sourceintercept, p_gamma_slope, NA),
           modelcomplexity = c(4,3,2,1))# assign a model complexity for sorting purposes

  all_aics <- bind_rows(gaussian_aics, gamma_aics)

  all_aics


}

#' Extract AIC for one model
#'
#' @param one_model a GLM from the list of models returned from fit_trend_models_biomass or fit_trend_models_energy
#'
#' @return data frame with model family, link function, formula, dataset name, AIC, AICc
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom MuMIn AICc
pull_model_aic <- function(one_model) {

  one_model_aic <- data.frame(model_AIC = AIC(one_model),
                              model_AICc = MuMIn::AICc(one_model)) %>%
    dplyr::mutate(model_family = one_model$family$family,
           model_link = one_model$family$link,
           model_formula = toString(one_model$formula[3]),
           matssname = one_model$data$matssname[1])
  one_model_aic
}

#' Get predicted trends from all models
#'
#' Get predicted trends from all models in a list
#'
#' @param some_models list of models produced using fit_trend_models_biomass or fit_trend_models_energy
#'
#' @return dataframe of predicted changes
#' @export
#'
#' @importFrom dplyr bind_rows
all_models_predicted_change <- function(some_models) {

  predicted_changes <- lapply(some_models, model_predicted_change)

  predicted_changes <- dplyr::bind_rows(predicted_changes)
  predicted_changes
}


#' Get predicted long-term change from a single model
#'
#' Primarily, calculates the long-term change (ratio of last fitted value to first fitted value) from the predicted values of a model.
#' Using this rather than the estimate of the slope terms, because the estimated slope terms from a Gamma model are not intuitive to interpret from the readout alone.
#'
#' @param one_model one model from the lists produced by fit_trend_models_biomass or fit_trend_models_energy
#'
#' @return data frame of predicted change, timespan, begin year, end year, correlation coefficient predicted/actual data, correlation coefficient squared, model info
#' @export
#'
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_wider
model_predicted_change <- function(one_model) {

  # Get predicted response values from the model
  model_dat <- one_model$data %>%
    mutate(predvals = predict(one_model, type = "response"))

  # Correlation coefficient to calculate a kind of r2
  model_cor <- cor(model_dat$predvals, one_model$model[,1])

  cor_squared = model_cor ^ 2


  # Get begin and end years and the time span
  model_begin <- min(model_dat$timeperiod)
  model_end <- max(model_dat$timeperiod)
  model_span <- model_end - model_begin


  model_change <- model_dat %>%
    dplyr::filter(timeperiod  %in% c(model_end, model_begin)) %>%
    dplyr::mutate(timeperiod_name = ifelse(timeperiod == model_end, "end", "begin")) %>%
    dplyr::select(timeperiod_name, predvals, source) %>%
    tidyr::pivot_wider(names_from = c(timeperiod_name, source), values_from = predvals) %>%
    dplyr::mutate(ratio_sim = end_sim / begin_sim,
           ratio_real = end_real / begin_real,
           model_cor = model_cor,
           model_cor_squared = cor_squared) %>%
    dplyr::mutate(model_family = one_model$family$family,
           model_link = one_model$family$link,
           model_formula = toString(one_model$formula[3]),
           model_begin = model_begin,
           model_end = model_end,
           model_span = model_span,
           matssname = model_dat$matssname[1],
           model_observations = length(unique(model_dat$timeperiod))
    )

  model_change
}





#' Compare sub-community properties
#'
#' Compares the mean body size, ISDs, and species composition for the first to the last 5 years of the timeseries.
#'
#' @param dat a MATSS dataset
#'
#' @return df
#' @export
#'
#' @importFrom dplyr filter left_join select mutate group_by_all ungroup group_by summarize
#' @importFrom tidyr pivot_wider
#' @importFrom vegan vegdist
compare_community_structure <- function(dat) {

  dat_focal_years <-pull_focal_years(dat)

  dat_null_mod <- relabund_null_model(dat_focal_years, resp_seed = 1989)

  dat_isd_actual <- simulate_isd_ts(dat_focal_years, isd_seed = 1989)

  dat_isd_null <- simulate_isd_ts(dat_null_mod, isd_seed = 1989)


  # pull first and last 5 years
  years <- unique(dat_focal_years$covariates$year)

  begin_years <- years[1:5]
  end_years <- years[(length(years) - 4):length(years)]

  # pull isds (mass estimates for each individual) for begin and end time periods
  begin_isd_real <- dplyr::filter(dat_isd_actual$isd, year %in% begin_years)
  end_isd_real <- dplyr::filter(dat_isd_actual$isd, year %in% end_years )


  # ISD overlap
  # Get density smooths by fitting a Gaussian mixture model to the ISD, extracting the density smooth, and scaling so the area under the density smooth sums to 1
  # Then combine into a df so we can compare begin to end
  begin_isd_gmm_real <- add_gmm(begin_isd_real)
  end_isd_gmm_real <- add_gmm(end_isd_real) %>%
    dplyr::rename(density_end = density)

  real_gmms <- dplyr::left_join(begin_isd_gmm_real, end_isd_gmm_real)

  # Calculate overlap metric between the two density smooths, and turnover as 1-overlap

  real_overlap <- real_gmms %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(minDensity = min(density, density_end)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(overlap = sum(minDensity)) %>%
    dplyr::mutate(turnover = 1 - overlap)


  # Repeat the above for the null model
  begin_isd_sim <- filter(dat_isd_null$isd, year %in% begin_years)
  end_isd_sim <- filter(dat_isd_null$isd, year %in% end_years)

  begin_isd_gmm_sim <- add_gmm(begin_isd_sim)
  end_isd_gmm_sim <- add_gmm(end_isd_sim) %>%
    rename(density_end = density)

  sim_gmms <- left_join(begin_isd_gmm_sim, end_isd_gmm_sim)

  sim_overlap <- sim_gmms %>%
    dplyr::group_by_all() %>%
   dplyr::mutate(minDensity = min(density, density_end)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(overlap = sum(minDensity))%>%
    dplyr::mutate(turnover = 1 - overlap)

  # Calcualte mean mass in first and last 5 year time periods for real and null model
  real_begin_mean_mass <- mean(begin_isd_real$mass)
  real_end_mean_mass <- mean(end_isd_real$mass)
  sim_begin_mean_mass <- mean(begin_isd_sim$mass)
  sim_end_mean_mass <- mean(end_isd_sim$mass)


  # For species compositional turnover (BCD)
  # Create species abundance matrix for first 5/last 5 years and calculate BCD

  sp_matrix <- bind_rows(begin_isd_real, end_isd_real) %>%
    dplyr::mutate(timeperiod = ifelse(year < 2000, "begin", 'end')) %>%
    dplyr::group_by(timeperiod, id) %>%
    dplyr::summarize(abund = dplyr::n()) %>%
    tidyr::pivot_wider(id_cols= timeperiod, names_from = id, values_from = abund, values_fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-timeperiod) %>%
    as.matrix()

  real_bcd <- vegan::vegdist(sp_matrix, "bray")[1]

  out <- data.frame(
    real_overlap = real_overlap$overlap[1],
    sim_overlap = sim_overlap$overlap[1],
    real_turnover = real_overlap$turnover[1],
    sim_turnover = sim_overlap$turnover[1],
    begin_years = toString(begin_years),
    end_years = toString(end_years),
    real_begin_mean_mass = real_begin_mean_mass,
    real_end_mean_mass = real_end_mean_mass,
    sim_begin_mean_mass = sim_begin_mean_mass,
    sim_end_mean_mass = sim_end_mean_mass,
    real_mass_ratio = real_end_mean_mass / real_begin_mean_mass,
    sim_mass_ratio = sim_end_mean_mass / sim_begin_mean_mass,
    sp_turnover_bcd = real_bcd
  ) %>%
    dplyr::bind_cols(dat$metadata$location) %>%
    dplyr::mutate(matssname = paste0("bbs_rtrg_", route, "_", statenum))

  out
}

#' Fit a GMM
#'
#' @param size_vect vector of sizes (logged or not)
#' @param max_G number of Gaussians
#'
#' @return GMM fit with up to max_G Gaussians
#' @export
#' @importFrom mclust densityMclust
fit_gmm <- function(size_vect, max_G = 15) {
  library(mclust)
  this_gmm <- mclust::densityMclust(size_vect, G = c(1:max_G) )

  return(this_gmm)

}

#' Add GMM density to an ISD
#'
#' Fits to log of mass
#'
#' @param isd from simulate_isd_ts
#' @param max_size default 15000 but changeable for mammals
#' @param max_G max n gaussians
#'
#' @return dataframe with columns mass, density
#' @export
#'
#' @importFrom dplyr mutate
add_gmm <- function(isd, max_size = 15000, max_G = 15) {

  isd <- isd %>%
    dplyr::mutate(logmass = log(mass))

  gmm <- fit_gmm(isd$logmass, max_G)


  gmm_isd <- data.frame(logmass = seq(0, log(max_size), length.out = 1000))
  gmm_isd$dens <- predict(gmm, newdata = gmm_isd$logmass)


  isd_gmm <- data.frame(
    mass = gmm_isd$logmass,
    density = (gmm_isd$dens)/ sum(gmm_isd$dens)
  )

  return(isd_gmm)

}


#' Pull years
#'
#' @param dataset a MATSS dataset
#'
#' @return df with year coverage
#' @export
#'

get_years <- function(dataset) {

  yeardat <- data.frame(
    year = dataset$covariates$year,
    route = dataset$metadata$route,
    region = dataset$metadata$region,
    location.bcr = dataset$metadata$location$bcr,
    location.statenum = dataset$metadata$location$statenum,
    location.routename = dataset$metadata$location$routename
  )

  return(yeardat)

}


