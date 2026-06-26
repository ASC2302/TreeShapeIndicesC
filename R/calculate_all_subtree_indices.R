#' Calculate tree shape indices for all valid subtrees
#'
#' This function reads a tree or collection of trees, extracts all subtrees
#' within a specified tip-count range, and calculates tree shape indices for
#' each subtree using [all_indices()].
#'
#' Invalid trees are rejected before calculation. Invalid outputs, such as J
#' values outside [0, 1], are also rejected.
#'
#' @param file A file path to a tree file in Newick or NEXUS format, a Newick
#'   string, or a phylo/multiPhylo object.
#' @param min_tips Minimum number of tips for subtrees to be included.
#' @param max_tips Maximum number of tips for subtrees to be included.
#' @param include_full_tree Logical. Should the full tree also be included in
#'   the results?
#'
#' @return A data frame containing calculated indices for each subtree, along
#'   with metadata such as tree number, subtree number, root label, leaf count,
#'   and calculation time.
#'
#' @export
calculate_all_subtree_indices <- function(
  file,
  min_tips = 0,
  max_tips = 20000,
  include_full_tree = FALSE
) {

  tree_obj <- read_convert(file)
  tree_list <- if (inherits(tree_obj, "multiPhylo")) tree_obj else list(tree_obj)

  results_list <- list()
  row_counter <- 1

  tictoc::tic("Total subtree processing time")

  for (tree_i in seq_along(tree_list)) {
    message("Processing full tree ", tree_i, " of ", length(tree_list))

    tree <- normalize_phylo(
      tree_list[[tree_i]],
      context = paste0("tree ", tree_i)
    )

    if (!inherits(tree, "phylo")) {
      warning(paste("Skipping tree", tree_i, "because it is not a valid phylo object."))
      next
    }

    if (ape::Ntip(tree) == 0) {
      warning(paste("Skipping tree", tree_i, "because it has no tips."))
      next
    }

    sts <- tryCatch(
      suppressWarnings(ape::subtrees(tree)),
      error = function(e) {
        warning(paste("Could not generate subtrees for tree", tree_i, ":", e$message))
        return(list())
      }
    )

    sts <- Filter(function(x) !is.null(x) && inherits(x, "phylo"), sts)

    sts <- lapply(
      seq_along(sts),
      function(i) {
        normalize_phylo(
          sts[[i]],
          context = paste0("tree ", tree_i, ", generated subtree ", i)
        )
      }
    )

    if (length(sts) > 0) {
      keep_idx <- vapply(
        sts,
        function(x) {
          nt <- ape::Ntip(x)
          nt >= min_tips && nt <= max_tips
        },
        logical(1)
      )

      sts <- sts[keep_idx]
    }

    if (include_full_tree) {
      sts <- c(list(tree), sts)
    }

    if (length(sts) > 0) {
      tree_strings <- vapply(sts, safe_write_tree, character(1))
      sts <- sts[!duplicated(tree_strings)]
    }

    message("  Number of trees/subtrees to process: ", length(sts))

    if (length(sts) == 0) next

    for (sub_i in seq_along(sts)) {
      st <- normalize_phylo(
        sts[[sub_i]],
        context = paste0("tree ", tree_i, ", subtree ", sub_i)
      )

      if (!inherits(st, "phylo")) next
      if (ape::Ntip(st) == 0) next

      timing <- system.time({
        idx <- tryCatch(
          {
            idx <- all_indices(st)

            validate_index_output(
              idx,
              context = paste0("tree ", tree_i, ", subtree ", sub_i)
            )

            idx
          },
          error = function(e) {
            warning(
              paste(
                "Error for tree",
                tree_i,
                "subtree",
                sub_i,
                ":",
                e$message
              )
            )
            return(NULL)
          }
        )
      })

      if (is.null(idx)) next

      root_label <- if (!is.null(st$node.label) && length(st$node.label) > 0) {
        as.character(st$node.label[1])
      } else {
        NA_character_
      }

      results_list[[row_counter]] <- data.frame(
        TreeNumber = tree_i,
        SubtreeNumber = sub_i,
        RootLabel = root_label,
        LeafCount = ape::Ntip(st),
        D0N = safe_extract(idx$D0N),
        D1N = safe_extract(idx$D1N),
        J1N = safe_extract(idx$J1N),
        D0S = safe_extract(idx$D0S),
        D1S = safe_extract(idx$D1S),
        J1S = safe_extract(idx$J1S),
        D0L = safe_extract(idx$D0L),
        D1L = safe_extract(idx$D1L),
        J1L = safe_extract(idx$J1L),
        TimeTakenSeconds = unname(timing["elapsed"]),
        stringsAsFactors = FALSE
      )

      row_counter <- row_counter + 1

      if (sub_i %% 25 == 0) {
        message("   Processed ", sub_i, " subtrees")
      }
    }
  }

  tictoc::toc()

  if (length(results_list) == 0) {
    return(data.frame(
      TreeNumber = integer(),
      SubtreeNumber = integer(),
      RootLabel = character(),
      LeafCount = integer(),
      D0N = numeric(),
      D1N = numeric(),
      J1N = numeric(),
      D0S = numeric(),
      D1S = numeric(),
      J1S = numeric(),
      D0L = numeric(),
      D1L = numeric(),
      J1L = numeric(),
      TimeTakenSeconds = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  dplyr::bind_rows(results_list)
}
