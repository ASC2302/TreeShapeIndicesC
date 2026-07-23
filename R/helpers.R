# Validate a TRUE/FALSE option used by the public wrappers.
validate_logical_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be either TRUE or FALSE.", call. = FALSE)
  }

  x
}


# Validate a phylo object before calculation
validate_phylo <- function(
  tr,
  allow_zero_lengths = TRUE,
  allow_missing_lengths = TRUE,
  ignore_branch_lengths = FALSE,
  context = "tree"
) {
  ignore_branch_lengths <- validate_logical_flag(
    ignore_branch_lengths,
    "ignore_branch_lengths"
  )

  if (!inherits(tr, "phylo")) {
    stop(context, ": input is not a valid phylo object.", call. = FALSE)
  }

  if (is.null(tr$edge) || nrow(tr$edge) == 0) {
    stop(
      context,
      ": tree has no edges. It may be empty or have fewer than two tips.",
      call. = FALSE
    )
  }

  n_tips <- tryCatch(ape::Ntip(tr), error = function(e) NA_integer_)

  if (is.na(n_tips) || n_tips < 2) {
    stop(context, ": tree must have at least two tips.", call. = FALSE)
  }

  if (is.null(tr$tip.label) || length(tr$tip.label) != n_tips) {
    stop(context, ": tree has missing or invalid tip labels.", call. = FALSE)
  }

  if (anyDuplicated(tr$tip.label)) {
    duplicated_tips <- unique(tr$tip.label[duplicated(tr$tip.label)])

    stop(
      context,
      ": tree contains duplicate tip labels: ",
      paste(duplicated_tips, collapse = ", "),
      call. = FALSE
    )
  }

  # The topology is still validated, but supplied edge-length values are not
  # inspected because normalize_phylo() will replace all of them with 1.
  if (ignore_branch_lengths) {
    return(TRUE)
  }

  if (is.null(tr$edge.length) || length(tr$edge.length) == 0) {
    if (allow_missing_lengths) {
      return(TRUE)
    }

    stop(context, ": tree is missing branch lengths.", call. = FALSE)
  }

  if (!is.numeric(tr$edge.length)) {
    stop(context, ": branch lengths must be numeric.", call. = FALSE)
  }

  if (length(tr$edge.length) != nrow(tr$edge)) {
    stop(
      context,
      ": number of branch lengths does not match number of edges. ",
      "Expected ", nrow(tr$edge), " but found ", length(tr$edge.length), ".",
      call. = FALSE
    )
  }

  edge_lengths <- tr$edge.length

  if (any(is.nan(edge_lengths))) {
    bad <- which(is.nan(edge_lengths))

    stop(
      context,
      ": tree contains NaN branch length(s) at edge position(s): ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " ...",
      call. = FALSE
    )
  }

  if (any(is.na(edge_lengths))) {
    bad <- which(is.na(edge_lengths))

    stop(
      context,
      ": tree contains NA branch length(s) at edge position(s): ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " ...",
      call. = FALSE
    )
  }

  if (any(is.infinite(edge_lengths))) {
    bad <- which(is.infinite(edge_lengths))

    stop(
      context,
      ": tree contains infinite branch length(s) at edge position(s): ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " ...",
      call. = FALSE
    )
  }

  if (any(edge_lengths < 0)) {
    bad <- which(edge_lengths < 0)

    stop(
      context,
      ": tree contains ",
      length(bad),
      " negative branch length(s). ",
      "Branch lengths must be non-negative. First bad edge position(s): ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " ...",
      call. = FALSE
    )
  }

  if (!allow_zero_lengths && any(edge_lengths == 0)) {
    bad <- which(edge_lengths == 0)

    stop(
      context,
      ": tree contains ",
      length(bad),
      " zero-length branch(es). ",
      "Zero branch lengths are not allowed under the current settings.",
      call. = FALSE
    )
  }

  TRUE
}


# Validate output from all_indices()
validate_index_output <- function(
  idx,
  context = "index output",
  tolerance = 1e-8
) {
  required_names <- c(
    "D0N", "D1N", "J1N",
    "D0S", "D1S", "J1S",
    "D0L", "D1L", "J1L"
  )

  missing_names <- setdiff(required_names, names(idx))

  if (length(missing_names) > 0) {
    stop(
      context,
      ": index result is missing required value(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  values <- vapply(
    required_names,
    function(nm) {
      x <- idx[[nm]]

      if (is.null(x) || length(x) != 1) {
        stop(
          context,
          ": ",
          nm,
          " must be a single numeric value.",
          call. = FALSE
        )
      }

      as.numeric(x)
    },
    numeric(1)
  )

  if (any(!is.finite(values))) {
    bad <- names(values)[!is.finite(values)]

    stop(
      context,
      ": non-finite index value(s) produced: ",
      paste(paste0(bad, "=", values[bad]), collapse = ", "),
      call. = FALSE
    )
  }

  j_names <- c("J1N", "J1S", "J1L")
  j_values <- values[j_names]

  bad_j <- j_names[j_values < -tolerance | j_values > 1 + tolerance]

  if (length(bad_j) > 0) {
    stop(
      context,
      ": J index value(s) outside expected range [0, 1]: ",
      paste(paste0(bad_j, "=", signif(j_values[bad_j], 10)), collapse = ", "),
      call. = FALSE
    )
  }

  d_names <- c("D0N", "D1N", "D0S", "D1S", "D0L", "D1L")
  d_values <- values[d_names]

  bad_d <- d_names[d_values <= 0]

  if (length(bad_d) > 0) {
    stop(
      context,
      ": D index value(s) must be positive: ",
      paste(paste0(bad_d, "=", signif(d_values[bad_d], 10)), collapse = ", "),
      call. = FALSE
    )
  }

  TRUE
}


# Make a phylo object safe to use/write:
# - reject nonsense branch lengths unless they are deliberately ignored
# - replace missing or ignored branch lengths with unit lengths
# - remove invalid node labels
normalize_phylo <- function(
  tr,
  allow_zero_lengths = TRUE,
  allow_missing_lengths = TRUE,
  ignore_branch_lengths = FALSE,
  context = "tree"
) {
  if (!inherits(tr, "phylo")) return(tr)

  ignore_branch_lengths <- validate_logical_flag(
    ignore_branch_lengths,
    "ignore_branch_lengths"
  )

  validate_phylo(
    tr,
    allow_zero_lengths = allow_zero_lengths,
    allow_missing_lengths = allow_missing_lengths,
    ignore_branch_lengths = ignore_branch_lengths,
    context = context
  )

  if (ignore_branch_lengths ||
      is.null(tr$edge.length) ||
      length(tr$edge.length) == 0) {
    tr$edge.length <- rep(1, nrow(tr$edge))
  }

  if (!is.null(tr$node.label) && length(tr$node.label) != tr$Nnode) {
    tr$node.label <- NULL
  }

  tr
}


resolve_node_abundances <- function(
  node_abundances = NULL,
  ignore_node_sizes = FALSE
) {
  ignore_node_sizes <- validate_logical_flag(
    ignore_node_sizes,
    "ignore_node_sizes"
  )

  if (ignore_node_sizes) NULL else node_abundances
}


prepare_index_inputs <- function(
  file,
  node_abundances = NULL,
  ignore_branch_lengths = FALSE,
  ignore_node_sizes = FALSE
) {
  ignore_branch_lengths <- validate_logical_flag(
    ignore_branch_lengths,
    "ignore_branch_lengths"
  )

  tree <- assert_single_phylo(
    read_convert(
      file,
      ignore_branch_lengths = ignore_branch_lengths
    )
  )

  list(
    tree = tree,
    node_abundances = resolve_node_abundances(
      node_abundances,
      ignore_node_sizes = ignore_node_sizes
    )
  )
}


safe_extract <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  as.numeric(x)
}


safe_write_tree <- function(tr, ignore_branch_lengths = FALSE) {
  tr <- normalize_phylo(
    tr,
    ignore_branch_lengths = ignore_branch_lengths
  )
  ape::write.tree(tr)
}


assert_single_phylo <- function(tree, arg = "file") {
  if (inherits(tree, "multiPhylo")) {
    stop(
      "`", arg, "` contains multiple trees. Use calculate_all_subtree_indices() or process_tree_folder() for multiPhylo inputs.",
      call. = FALSE
    )
  }

  if (!inherits(tree, "phylo")) {
    stop("`", arg, "` must be a phylo object after conversion.", call. = FALSE)
  }

  tree
}
