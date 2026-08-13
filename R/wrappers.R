#' Calculate all universal tree shape indices for one tree
#'
#' This wrapper validates the input tree and checks that the output indices are
#' within expected ranges. In particular, J indices are checked to be in [0, 1].
#'
#' @param file A file path, Newick string, or phylo object.
#' @param node_abundances Optional node abundance data frame passed to the C++
#'   index functions.
#' @param validate_output Logical. Should output indices be checked?
#' @param ignore_branch_lengths Logical. Replace all edge lengths with 1 before
#'   calculation, producing a topology/unit-edge result.
#' @param ignore_node_sizes Logical. Ignore `node_abundances` and use the C++
#'   default of equal abundance on tips and zero direct abundance on internal
#'   nodes.
#'
#' @return A list containing D0N, D1N, J1N, D0S, D1S, J1S, D0L, D1L, and J1L.
#'
#' @examples
#' tree_file <- tree_index_example("rphylo.nwk")
#' tree <- read_convert(tree_file)
#'
#' # Use a small section of the bundled tree so the example runs quickly.
#' example_tree <- ape::keep.tip(tree, tree$tip.label[seq_len(10)])
#'
#' result <- all_indices(
#'   example_tree,
#'   ignore_branch_lengths = TRUE
#' )
#'
#' result$J1N
#'
#' @export
all_indices <- function(
  file,
  node_abundances = NULL,
  validate_output = TRUE,
  ignore_branch_lengths = FALSE,
  ignore_node_sizes = FALSE
) {
  inputs <- prepare_index_inputs(
    file = file,
    node_abundances = node_abundances,
    ignore_branch_lengths = ignore_branch_lengths,
    ignore_node_sizes = ignore_node_sizes
  )

  idx <- all_indices_cpp(inputs$tree, inputs$node_abundances)

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
#' @param ignore_branch_lengths Logical. Replace all edge lengths with 1 before
#'   calculation.
#' @param ignore_node_sizes Logical. Ignore `node_abundances` and use equal tip
#'   abundance with zero direct internal-node abundance.
#'
#' @return A list or numeric value depending on `individual`.
#'
#' @export
node <- function(
  file,
  node_abundances = NULL,
  index_letter = "D",
  q = 1,
  individual = FALSE,
  ignore_branch_lengths = FALSE,
  ignore_node_sizes = FALSE
) {
  inputs <- prepare_index_inputs(
    file = file,
    node_abundances = node_abundances,
    ignore_branch_lengths = ignore_branch_lengths,
    ignore_node_sizes = ignore_node_sizes
  )

  node_cpp(
    inputs$tree,
    inputs$node_abundances,
    index_letter,
    q,
    individual
  )
}

#' Calculate star or longitudinal tree shape indices for one tree
#'
#' @param file A file path, Newick string, or phylo object.
#' @param node_abundances Optional node abundance data frame.
#' @param mean_type Either "Star" or "Longitudinal".
#' @param index_letter Either "D" or "J" when `individual = TRUE`.
#' @param q Diversity order used when `individual = TRUE`.
#' @param individual Logical. Return one individual index instead of all indices.
#' @param ignore_branch_lengths Logical. Replace all edge lengths with 1 before
#'   calculation.
#' @param ignore_node_sizes Logical. Ignore `node_abundances` and use equal tip
#'   abundance with zero direct internal-node abundance.
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
  individual = FALSE,
  ignore_branch_lengths = FALSE,
  ignore_node_sizes = FALSE
) {
  inputs <- prepare_index_inputs(
    file = file,
    node_abundances = node_abundances,
    ignore_branch_lengths = ignore_branch_lengths,
    ignore_node_sizes = ignore_node_sizes
  )

  long_star_cpp(
    inputs$tree,
    inputs$node_abundances,
    mean_type,
    index_letter,
    q,
    individual
  )
}
