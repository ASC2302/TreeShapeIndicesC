#' Read and convert a tree file into a phylo or multiPhylo object
#'
#' This function attempts to read a tree from a file path, trying Newick and
#' NEXUS formats. If the input is already a phylo or multiPhylo object, it is
#' returned as is. The function also ensures that edge lengths are present,
#' assigning a default length of 1 if they are missing.
#'
#' @param file A file path to a tree file in Newick or NEXUS format, a Newick
#'   string, or a phylo/multiPhylo object.
#'
#' @return A phylo or multiPhylo object.
#'
#' @export
read_convert <- function(file) {

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

  input_context <- if (is.character(file) && length(file) == 1 && file.exists(file)) {
    paste0("input tree: ", basename(file))
  } else {
    "input tree"
  }

  if (inherits(tree, "phylo")) {
    tree <- normalize_phylo(
      tree,
      context = input_context
    )

  } else if (inherits(tree, "multiPhylo")) {
    tree <- lapply(
      seq_along(tree),
      function(i) {
        normalize_phylo(
          tree[[i]],
          context = paste0(input_context, ", tree ", i)
        )
      }
    )

    class(tree) <- "multiPhylo"
  }

  tree
}
