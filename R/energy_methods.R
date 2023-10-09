#' Estimate b (energy)
#'
#' @param body_mass mass
#' @return estimated b (metabolic rate) using pars from Fristoe 2015
#' @export

estimate_b <- function(body_mass) {
  return(10.5 * (body_mass ^ .713))
}

