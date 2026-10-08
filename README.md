# TreeShapeIndicesC

`TreeShapeIndicesC` calculates universal tree shape indices for phylogenetic
trees and subtrees using R and Rcpp.

## Installation from GitHub

```r
install.packages("remotes")
remotes::install_github("ASC2302/RUIindices2")
```

## Bundled example data

The package includes two raw Newick files under `inst/extdata/`:

- `rphylo.nwk`, generated using `ape::rphylo(500, 1, 0)`;
- `sim-taxa.nwk`, generated using
  `TreeSimGM::sim.taxa(1, 500, waitsp = "rexp(1.2)")[[1]]`.

Each file contains one rooted, fully bifurcating tree with 500 tips, unique tip
labels, and positive branch lengths. Use `tree_index_example()` to locate the
installed files:

```r
library(TreeShapeIndicesC)

tree_index_example()

rphylo_file <- tree_index_example("rphylo.nwk")
rphylo_tree <- read_convert(rphylo_file)
ape::Ntip(rphylo_tree)
```

A quick calculation can be demonstrated on a small section of either tree:

```r
example_tree <- ape::keep.tip(
  rphylo_tree,
  rphylo_tree$tip.label[seq_len(10)]
)

idx <- all_indices(
  example_tree,
  ignore_branch_lengths = TRUE
)

idx$J1N
```

The original files remain available so that the package's tree-reading
interface can be tested as well as its calculations.

## Ignoring branch lengths and node sizes

The public calculation functions accept two independent logical options:

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
