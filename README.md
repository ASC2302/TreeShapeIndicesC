# TreeShapeIndices

`TreeShapeIndices` calculates universal tree shape indices for phylogenetic trees and subtrees using R and Rcpp.

## Installation from GitHub

```r
install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USERNAME/TreeShapeIndices")
```

Replace `YOUR-GITHUB-USERNAME` with the GitHub account or organisation that hosts the repository.

## Basic usage

```r
library(TreeShapeIndices)

results <- process_tree_folder(
  tree_folder = "C:/path/to/tree/files",
  min_tips = 0,
  max_tips = 20000,
  include_full_tree = FALSE
)
```

This writes one CSV file per input tree file. By default, files ending in `.txt`, `.nwk`, `.tree`, `.tre`, `.nex`, or `.nexus` are processed.

## Process one tree

```r
results <- calculate_all_subtree_indices(
  file = "C:/path/to/tree.nwk",
  min_tips = 100,
  max_tips = 1000,
  include_full_tree = FALSE
)
```

## Calculate indices for one full tree

```r
idx <- all_indices("C:/path/to/tree.nwk")
idx
```

## Development

After changing the C++ exports, regenerate the Rcpp glue files with:

```r
Rcpp::compileAttributes()
devtools::document()
```

Then test installation locally with:

```r
devtools::check()
```

## C++ source

The main Rcpp source file is stored at:

```text
src/main_pipeline.cpp
```

After editing exported C++ functions, regenerate the Rcpp export files with:

```r
Rcpp::compileAttributes()
devtools::document()
```

