# Applying the complete TreeShapeIndicesC patch

This archive is an overlay for the `TreeShapeIndicesC` package. It includes the
branch-length/node-size options, bundled example data, tests, and examples.

## Copy the files

From the extracted archive, copy the following into the package repository,
allowing matching files to be replaced:

- `R/` to the package's `R/` directory;
- `src/RcppExports.cpp` to `src/RcppExports.cpp`;
- `inst/` to the package root;
- `tests/` to the package root;
- `DESCRIPTION`, `NAMESPACE`, `README.md`, and `PATCH_NOTES.md` to the package
  root.

The raw data must finish at:

```text
inst/extdata/rphylo.nwk
inst/extdata/sim-taxa.nwk
```

Do not copy these files into `data/`; they are raw Newick inputs and belong in
`inst/extdata/`.

## Regenerate and check the package

Open the package project in RStudio and run from the package root:

```r
Rcpp::compileAttributes()
devtools::document()
devtools::test()
devtools::check()
devtools::install()
```

`Rcpp::compileAttributes()` is important because older generated files in the
repository used the `_TreeShapeIndices_` prefix instead of
`_TreeShapeIndicesC_`.

`devtools::document()` regenerates the `NAMESPACE` and help files, including the
export and documentation for `tree_index_example()`.

## Confirm the example files after installation

```r
library(TreeShapeIndicesC)

tree_index_example()
file.exists(tree_index_example("rphylo.nwk"))
file.exists(tree_index_example("sim-taxa.nwk"))
```

The expected result is two `TRUE` values.
