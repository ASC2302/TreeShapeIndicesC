#' Read and convert a tree file into a phylo or multiPhylo object
#'
#' This function attempts to read a tree from a file path, trying Newick and
#' NEXUS formats. If the input is already a phylo or multiPhylo object, it is
#' returned after validation. Missing branch lengths are assigned a default
#' value of 1. When `ignore_branch_lengths = TRUE`, all supplied branch lengths
#' are replaced with 1 before the tree is used by the index calculations.
#'
#' @param file A file path to a tree file in Newick or NEXUS format, a Newick
#'   string, or a phylo/multiPhylo object.
#' @param ignore_branch_lengths Logical. Replace every edge length with 1 and
#'   calculate from topology rather than the supplied branch lengths.
#'
#' @return A phylo or multiPhylo object.
#'
#' @export
read_convert <- function(file, ignore_branch_lengths = FALSE) {
  ignore_branch_lengths <- validate_logical_flag(
    ignore_branch_lengths,
    "ignore_branch_lengths"
  )

  if (inherits(file, "phylo") || inherits(file, "multiPhylo")) {
    tree <- file
  } else {
    suppressWarnings({
      tree <- try(ape::read.tree(file), silent = TRUE)

      if (inherits(tree, "try-error")) {
        tree <- try(ape::read.nexus(file), silent = TRUE)
      }

      if (inherits(tree, "try-error")) {
        tree <- try(ape::read.tree(text = file), silent = TRUE)
      }
    })

    if (inherits(tree, "try-error")) {
      stop(
        "Tree must be in Newick or NEXUS format, ",
        "or be a phylo/multiPhylo object.",
        call. = FALSE
      )
    }
  }

  input_context <- if (
    is.character(file) && length(file) == 1 && file.exists(file)
  ) {
    paste0("input tree: ", basename(file))
  } else {
    "input tree"
  }

  if (inherits(tree, "phylo")) {
    tree <- normalize_phylo(
      tree,
      ignore_branch_lengths = ignore_branch_lengths,
      context = input_context
    )
  } else if (inherits(tree, "multiPhylo")) {
    tree <- lapply(
      seq_along(tree),
      function(i) {
        normalize_phylo(
          tree[[i]],
          ignore_branch_lengths = ignore_branch_lengths,
          context = paste0(input_context, ", tree ", i)
        )
      }
    )

    class(tree) <- "multiPhylo"
  }

  tree
}
