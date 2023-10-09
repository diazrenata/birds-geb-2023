#' Records of mean and standard deviation masses for species in BBS, from Dunning (2008)
#'
#' From Dunning, raw records of mean and (where present) standard deviation for species in the BBS.
#'
#' @format A data frame with 897 rows and 19 variables:
#' \describe{
#'   \item{species_id}{also called aou, species code for BBS}
#'   \item{id}{species ID column}
#'   \item{english_common_name}{inherited}
#'   \item{sporder}{inherited}
#'   \item{family}{inherited}
#'   \item{genus}{Genus}
#'   \item{species}{Species}
#'   \item{mass}{mass (in g) from Dunning (2008)}
#'   \item{sd}{sd of mass from Dunning, if reported}
#'   \item{sex}{sex for record from Dunning, if reported}
#'   \item{subspecies}{subspecies from Dunning, if reported}
#'   \item{location}{location from Dunning, if reported}
#'   \item{name_mismatch}{1 or NA, whether the name in Dunning does not match the name in BBS. Typically slight spelling differences.}
#'   \item{name_notes}{name as listed in Dunning, if different from that listed in BBS}
#'   \item{not_in_dunning}{1 or NA, whether this species is not found in Dunning}
#'   \item{close_genus}{if a species is not in Dunning (2008) because of a name change or species split since 2008, genus of new classification}
#'   \item{close_species}{new species classification}
#'   \item{close_subspecies}{new subspecies classification}
#'   \item{close_species_notes}{description of species reclassification, if original as named is not in Dunning}
#' }
"raw_masses"
