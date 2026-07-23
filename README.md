# TreeShapeIndicesC

`TreeShapeIndicesC` calculates universal tree shape indices for phylogenetic
trees and subtrees using R and Rcpp.

## Installation from GitHub

```r
install.packages("remotes")
remotes::install_github("ASC2302/TreeShapeIndicesC")
```

## Ignoring branch lengths and node sizes

The public calculation functions now accept two independent logical options:

```r
ignore_branch_lengths = FALSE
ignore_node_sizes = FALSE
```

- `ignore_branch_lengths = TRUE` replaces every edge length with `1`, so the
  result depends on topology/unit edges rather than the supplied lengths.
- `ignore_node_sizes = TRUE` discards `node_abundances`. The C++ default is
  then used: equal direct abundance on tips and zero direct abundance on
  internal nodes.

### Full tree

```r
library(TreeShapeIndicesC)

idx <- all_indices(
  "C:/path/to/tree.nwk",
  node_abundances = node_sizes,
  ignore_branch_lengths = TRUE,
  ignore_node_sizes = TRUE
)
```

### Subtrees

```r
results <- calculate_all_subtree_indices(
  file = "C:/path/to/tree.nwk",
  min_tips = 100,
  max_tips = 1000,
  include_full_tree = FALSE,
  node_abundances = node_sizes,
  ignore_branch_lengths = TRUE,
  ignore_node_sizes = FALSE
)
```

### Folder processing

```r
results <- process_tree_folder(
  tree_folder = "C:/path/to/tree/files",
  min_tips = 0,
  max_tips = 20000,
  include_full_tree = FALSE,
  ignore_branch_lengths = TRUE,
  ignore_node_sizes = TRUE
)
```

By default, files ending in `.txt`, `.nwk`, `.tree`, `.tre`, `.nex`, or
`.nexus` are processed.
