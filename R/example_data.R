#' Dataset for use in vignettes: example_data
#'
#' A dataset containing metatranscriptome data for responder taxon *Pseudoalteromonas sp.* and companion taxon *Chaetoceors dichaeta*.
#'
#' @format A data frame with 22 rows and 9 variables:
#' \describe{
#'   \item{sample_name_r}{Sample id with replicate information (character)}
#'   \item{xx}{Observed abundance of responder taxon (integer)}
#'   \item{xstar}{Observed abundance of companion taxon (integer)}
#'   \item{sample_name}{Sample id (character)}
#'   \item{temp}{Temperature in Celsius of ocean from which sample was taken (integer)}
#'   \item{K02520_counts_sum}{Observed abundance of gene K02520 expressed by the responder taxon (integer)}
#'   \item{K01006_counts_sum}{Observed abundance of gene K01006 expressed by the responder taxon (integer)}
#'   \item{K01497_counts_sum}{Observed abundance of gene K01497 expressed by the responder taxon (integer)}
#'   \item{K03106_counts_sum}{Observed abundance of gene K03106 expressed by the responder taxon (integer)}
#' }
#' @source Bartolek et al. "Functional patterns of microbial interaction in the North Pacific revealed through statistical modeling of environmental metatranscriptomes." (2025+)
"example_data"
