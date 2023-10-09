#' Simulate an ISD for a single year
#'
#' Just using the abundances.
#'
#' @param abundances long dataframe of `id`, `abundance`, and any other identifying columns (such as `year`)
#' @param sd_dat optionally pass table with `id`, `mean_mass`, `mean_sd`
#' @param use_th_pars default FALSE, use Thibault (2010) parameter estiamtes?
#' @param isd_seed optionally specify seed to use
#'
#' @return dataframe with one row per simulated individual, columns `id`, `mass`, `isd_seed`, and any columns in `abundances`
#' @export
#'
#' @importFrom dplyr select distinct mutate left_join bind_rows
simulate_isd <- function(abundances, sd_dat = NULL, use_th_pars = F, isd_seed = NULL) {

  if(is.null(sd_dat)) {
    load(here::here("data", "sd_table.rda"))
    sd_dat = sd_table
  }

  if(is.null(isd_seed)) {
    isd_seed = sample.int(100000000, size = 1)
  }

  desc_data <- abundances %>%
    dplyr::select(-id, -abundance) %>%
    dplyr::distinct() %>%
    dplyr::mutate(isd_seed = isd_seed)

  abundances <- dplyr::left_join(abundances, sd_dat)

  if(use_th_pars) {

    abundances <- abundances %>%
      dplyr::mutate(mean_sd = mean_sd_th)

  }

  species_sims <- list()

  set.seed(isd_seed)

  for(i in 1:nrow(abundances)) {
    if(abundances$abundance[i] > 0 ) {
      species_sims[[i]] <- data.frame(
        id = abundances$id[i],
        mass = rnorm(n = abundances$abundance[i], mean = abundances$mean_mass[i], sd = abundances$mean_sd[i]),
        stringsAsFactors = F)

      while(any(species_sims[[i]]$mass < 0)) {

        negative_sizes <- which(species_sims[[i]]$mass < 0)

        species_sims[[i]]$mass[negative_sizes] <-  rnorm(n = length(negative_sizes), mean = abundances$mean_mass[i], sd = abundances$mean_sd[i])

      }

    }
  }

  simulated_individuals <- dplyr::bind_rows(species_sims) %>%
    cbind(desc_data)

  return(simulated_individuals)
}

#' Simulate ISD data for a timeseries
#'
#' For one year or all years.
#'
#' @param dataset list of `abundance`, `covariates`, `metadata`. See `MATSS`.
#' @param sd_dat optional, provide a dataframe with `id`, `mean_mass`, `mean_sd`
#' @param censusyears optional, specify the year to simulate for. If not provided simulates for all years.
#' @param use_th_pars default F, use Thibault (2010) scaling parameters or re-estimated ones?
#' @param isd_seed optionally specify the seed to use when simulating the ISDs. If not provided one will be chosen and returned with output.
#'
#' @return list of `isd`, `covariates`, `metadata`
#' @export
#'
#' @importFrom dplyr mutate filter bind_rows
#' @importFrom tidyr pivot_longer
simulate_isd_ts <- function(dataset, sd_dat = NULL, censusyears = NULL, use_th_pars = F, isd_seed = NULL) {

  if(is.null(sd_dat)) {
    load(here::here("data", "sd_table.rda"))
    sd_dat = sd_table
  }


  abundances_long <- dataset$abundance %>%
    dplyr::mutate(year = dataset$covariates$year) %>%
    tidyr::pivot_longer(-year, names_to = "id", values_to = "abundance")

  if(!is.null(censusyears)) {
    abundances_long <- dplyr::filter(abundances_long, year %in% censusyears)
  }

  annual_isds <- list()

  initial_isd_seed = isd_seed

  for(i in 1:length(unique(abundances_long$year))) {

    thisyear_abunds <- dplyr::filter(abundances_long,
                                     year == unique(abundances_long$year)[i])

    if(sum(thisyear_abunds$abundance) == 0) {

      annual_isds[[i]] <- NA
    } else {

      annual_isds[[i]] <- simulate_isd(abundances = thisyear_abunds,
                                       sd_dat = sd_dat,
                                       use_th_pars = use_th_pars,
                                       isd_seed = initial_isd_seed + i)
    }
  }

  annual_isds <- (annual_isds)[ which(!is.na(annual_isds))]
  all_annual_isds <- dplyr::bind_rows(annual_isds) %>%
    dplyr::mutate(isd_seed = initial_isd_seed)

  return(list(
    isd = all_annual_isds,
    covariates = dataset$covariates,
    metadata = dataset$metadata
  ))
}
