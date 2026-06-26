#' Calculate all universal tree shape indices for one tree
#'
#' This wrapper validates the input tree and checks that the output indices are
#' within expected ranges. In particular, J indices are checked to be in [0, 1].
#'
#' @param file A file path, Newick string, or phylo object.
#' @param node_abundances Optional node abundance data frame passed to the C++
#'   index functions.
#' @param validate_output Logical. Should output indices be checked?
#'
#' @return A list containing D0N, D1N, J1N, D0S, D1S, J1S, D0L, D1L, and J1L.
#'
#' @export
all_indices <- function(
  file,
  node_abundances = NULL,
  validate_output = TRUE
) {
  tree <- assert_single_phylo(read_convert(file))
  idx <- all_indices_cpp(tree, node_abundances)

  if (validate_output) {
    validate_index_output(idx, context = "all_indices()")
  }

  idx
}

#' Calculate node-based tree shape indices for one tree
#'
#' @param file A file path, Newick string, or phylo object.
#' @param node_abundances Optional node abundance data frame.
#' @param index_letter Either "D" or "J" when `individual = TRUE`.
#' @param q Diversity order used when `individual = TRUE`.
#' @param individual Logical. Return one individual index instead of all node indices.
#'
#' @return A list or numeric value depending on `individual`.
#'
#' @export
node <- function(
  file,
  node_abundances = NULL,
  index_letter = "D",
  q = 1,
  individual = FALSE
) {
  tree <- assert_single_phylo(read_convert(file))
  node_cpp(tree, node_abundances, index_letter, q, individual)
}

#' Calculate star or longitudinal tree shape indices for one tree
#'
#' @param file A file path, Newick string, or phylo object.
#' @param node_abundances Optional node abundance data frame.
#' @param mean_type Either "Star" or "Longitudinal".
#' @param index_letter Either "D" or "J" when `individual = TRUE`.
#' @param q Diversity order used when `individual = TRUE`.
#' @param individual Logical. Return one individual index instead of all indices.
#'
#' @return A list or numeric value depending on `individual`.
#'
#' @export
long_star <- function(
  file,
  node_abundances = NULL,
  mean_type = "Star",
  index_letter = "D",
  q = 1,
  individual = FALSE
) {
  tree <- assert_single_phylo(read_convert(file))
  long_star_cpp(tree, node_abundances, mean_type, index_letter, q, individual)
}
