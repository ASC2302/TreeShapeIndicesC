# Applying the TreeShapeIndicesC input-options patch

1. Copy the contents of `modified/R` into the package's `R` directory.
2. Copy `modified/src/RcppExports.cpp` into the package's `src` directory.
3. Copy `modified/DESCRIPTION` and `modified/README.md` to the package root.
4. Optionally copy the `tests` directory.

From RStudio, set the working directory to the package root and run:

```r
Rcpp::compileAttributes()
devtools::document()
devtools::test()
devtools::install()
```

`Rcpp::compileAttributes()` is important for the current repository because its
checked-in generated symbols still use the old package prefix
`_TreeShapeIndices_`. Regenerating them under package name `TreeShapeIndicesC`
creates the `_TreeShapeIndicesC_` symbols expected by the installed DLL.

## New arguments

The following public functions accept both switches:

- `all_indices()`
- `node()`
- `long_star()`
- `calculate_all_subtree_indices()`
- `process_tree_folder()`

```r
ignore_branch_lengths = FALSE
ignore_node_sizes = FALSE
```

`read_convert()` also accepts `ignore_branch_lengths`.
