#' Get SD parameters from linear model
#'
#' @param raw_size_data raw size data
#'
#' @return list of slope and intercept from lm
#' @export
#'
#' @importFrom dplyr filter mutate
get_sd_parameters <- function(raw_size_data) {

  sp_for_sd <- dplyr::filter(raw_size_data,
                             !is.na(sd)) %>%
    dplyr::mutate(mass = as.numeric(mass),
                  sd = as.numeric(sd)) %>%
    dplyr::mutate(var = sd^2) %>%
    dplyr::mutate(log_m = log(mass),
                  log_var = log(var))

  sd_fit <- lm(sp_for_sd, formula = log_var ~ log_m)

  intercept <- exp(sd_fit$coefficients[[1]])
  slope <- sd_fit$coefficient[[2]]

  return(list(intercept = intercept,
              slope = slope))

}


#' Estimate SD given mean and optional pars
#'
#' @param sp_mean mean to estimate from
#' @param pars if NA, uses Thibault fit; if supplied, uses given slope and intercept
#'
#' @return estimated sd
#' @export
#'
estimate_sd <- function(sp_mean, pars = NA) {


 # if(is.list(pars)) {

    fitted_sd = sqrt(pars$intercept * (sp_mean ^ pars$slope))

  # } else {
  #
  #   fitted_sd = sqrt(.0055 * (sp_mean ^ 1.98))
  # }

  return(fitted_sd)

}

#' Clean raw size data
#'
#' Updates names of species whose names have changed since 2010.
#'
#' @param raw_size_data raw size data
#'
#' @return cleaned data
#' @export
#'
#' @importFrom dplyr select mutate filter bind_rows
clean_sp_size_data <- function(raw_size_data) {

  colnames(raw_size_data)[1] <- "species_id"

  sp_clean <- raw_size_data %>%
    dplyr::select(-english_common_name, -sporder, -family) %>%
    dplyr::mutate(mass = as.numeric(mass))

  name_change <- dplyr::filter(sp_clean, not_in_dunning == 1)

  sp_clean <- dplyr::filter(sp_clean, is.na(not_in_dunning)) %>%
    dplyr::mutate(added_flag = NA)

  for(i in 1:nrow(name_change)) {

    if(!is.na(name_change$close_subspecies[i])) {
      matched_rows <- dplyr::filter(sp_clean,
                                    genus == name_change$close_genus[i],
                                    species == name_change$close_species[i],
                                    subspecies == name_change$close_subspecies[i])
    } else {
      matched_rows <- dplyr::filter(sp_clean,
                                    genus == name_change$close_genus[i],
                                    species == name_change$close_species[i])
    }

    sp_to_add <- matched_rows %>%
      dplyr:: mutate(species_id = name_change$species_id[i],
                     id = name_change$id[i],
                     added_flag = 1)

    sp_clean <- dplyr::bind_rows(sp_clean, sp_to_add)

  }

  return(sp_clean)
}


#' Add estimates for SDs
#'
#' @param clean_size_data cleaned data
#' @param sd_pars parameters as list of slope, intercept
#'
#' @return clean_size_data with filled in sds + column to flag if an estimate was added
#' @export
#'
add_estimated_sds <- function(clean_size_data, sd_pars) {

  clean_size_data$estimated_sd <- FALSE

  for(i in 1:nrow(clean_size_data)) {

    if(is.na(clean_size_data$sd[i])) {
      clean_size_data$estimated_sd[i] =  TRUE
      clean_size_data$sd[i] = estimate_sd(clean_size_data$mass[i], pars = sd_pars)
    } else {
      clean_size_data$estimated_sd[i] = FALSE
    }

  }

  return(clean_size_data)

}

#' Get species' mean size and sd
#'
#' @param sd_dat Dat with estimated sds
#'
#' @return Summarized to species mean size and sd
#' @export
#'
#' @importFrom dplyr group_by summarize ungroup
get_sp_mean_size <- function(sd_dat) {

  sp_means <- sd_dat %>%
    dplyr::group_by(species_id, id, genus, species) %>%
    dplyr::summarize(mean_mass = mean(mass),
              mean_sd = mean(sd, na.rm = F),
              contains_estimates = any(estimated_sd)) %>%
    dplyr::ungroup()


}


#' Generate table of species-level means for the mean and standard deviation of body mass for species in the Breeding Bird Survey
#'
#' Goes from the `raw_masses` dataframe (included in `bbssize`) of records of species' mean and (where provided) standard deviation of body mass from the CRC Handbook (Dunning 2008) to a table of species-level means for the mean and standard deviation of body mass, incorporating estimates for missing standard deviation records and resolving taxonomic updates between the publication of the CRC Handbook and present releases of the Breeding Bird Survey dataset (Paradieck et al. 2019).
#'
#' @param raw_size_data the `raw_masses` dataframe
#'
#' @references \itemize{
#'  \item{Dunning, J. B. (2008). CRC handbook of avian body masses (2nd ed.). CRC Press.}
#'   \item{Pardieck, K. L., Ziolkowski, D. J., Lutmerding, M., Aponte, V., & Hudson, M.-A. (2019). North American Breeding Bird Survey Dataset 1966—2018, version 2018.0. U.S. Geological Survey. https://doi.org/10.5066/P9HE8XYJ}
#' }
#'
#' @return a dataframe of species-level means for mean body size and standard deviation of body size
#' @keywords internal
#'
generate_sd_table <- function(raw_size_data) {

  # Calculate scaling parameters
  fitted_pars <- get_sd_parameters(raw_size_data)

  # Resolve name mismatches
  clean_size_dat <- clean_sp_size_data(raw_size_data)

  # Add estimates for missing standard deviation records
  sd_size_dat <-add_estimated_sds(
    clean_size_data = clean_size_dat,
    sd_pars = fitted_pars
  )

  # Summarize to species-level means for the mean and standard deviation of body mass
  sp_mean_size_dat <- get_sp_mean_size(sd_size_dat)

  sp_mean_size_dat
}
