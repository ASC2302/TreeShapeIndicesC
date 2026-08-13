#' Locate bundled example tree files
#'
#' `TreeShapeIndicesC` includes two raw Newick files that can be used in
#' examples, tests, and exploratory analyses. This helper returns their names
#' or the complete path to one installed file.
#'
#' @param file `NULL`, `"rphylo.nwk"`, or `"sim-taxa.nwk"`. When `NULL`, the
#'   available filenames are returned. Otherwise, the complete installed path
#'   to the selected file is returned.
#'
#' @return A character vector of available filenames when `file = NULL`, or a
#'   single complete file path when a filename is supplied.
#'
#' @details
#' `rphylo.nwk` contains a 500-tip tree generated using
#' `ape::rphylo(500, 1, 0)`. `sim-taxa.nwk` contains a 500-tip tree generated
#' using `TreeSimGM::sim.taxa(1, 500, waitsp = "rexp(1.2)")[[1]]`.
#'
#' @examples
#' tree_index_example()
#'
#' tree_file <- tree_index_example("rphylo.nwk")
#' tree <- read_convert(tree_file)
#' ape::Ntip(tree)
#'
#' @export
#'
#' @references
#' Wickham, H. and Bryan, J. (2023). R Packages (2e), Chapter 7: Data.
#'
#' @seealso [read_convert()], [all_indices()]
tree_index_example <- function(file = NULL) {
  available_files <- c("rphylo.nwk", "sim-taxa.nwk")

  if (is.null(file)) {
    return(available_files)
  }

  if (
    !is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      !nzchar(file) ||
      !file %in% available_files
  ) {
    stop(
      "`file` must be NULL or one of: ",
      paste(available_files, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  system.file(
    "extdata",
    file,
    package = "TreeShapeIndicesC",
    mustWork = TRUE
  )
}
