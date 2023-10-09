#' Means and (estimated) standard deviations of body size for species in the BBS.
#'
#' A table of the means and estimated standard deviations of body size for species in the BBS dataset based on values in Dunning (2008). Standard deviations, where missing, are estimated based on scaling parameters based on the species present in this dataset. Species-level means are reported for both the mean mean mass and mean standard deviation.
#'
#' @format A data frame with 421 rows and 7 variables:
#' \describe{
#'   \item{species_id}{also called aou, species code for BBS}
#'   \item{id}{species ID column}
#'   \item{genus}{Genus}
#'   \item{species}{Species}
#'   \item{mean_mass}{Mean mass across all records for that species in Dunning (2008)}
#'   \item{mean_sd}{Mean standard deviation across all records from Dunning (2008). For records missing a standard deviation, estimated based on the scaling relationship derived from this dataset. Estimated standard deviations are calculated prior to taking species-level means.}
#'   \item{contains_estimates}{TRUE/FALSE whether or not the mean_sd column contains estimated standard deviations}
#' }
"sd_table"
