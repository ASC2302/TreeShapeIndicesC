# Patch notes

## Added public options

```r
ignore_branch_lengths = FALSE
ignore_node_sizes = FALSE
```

These are available in `all_indices()`, `node()`, `long_star()`,
`calculate_all_subtree_indices()`, and `process_tree_folder()`.
`read_convert()` accepts `ignore_branch_lengths`.

## Behaviour

- Ignored branch lengths are replaced with unit edge lengths (`1`). This avoids
  zero-length integrals while removing the influence of the supplied lengths.
- Ignored node sizes cause `node_abundances` to be passed to C++ as `NULL`.
  The existing C++ default then assigns equal direct abundance to tips and zero
  direct abundance to internal nodes.
- Both options default to `FALSE`, preserving existing calls.

## Additional repair

The repository's generated Rcpp files use the old `_TreeShapeIndices_` symbol
prefix although the package is named `TreeShapeIndicesC`. The included generated
files use `_TreeShapeIndicesC_`. Running `Rcpp::compileAttributes()` from the
package root will regenerate the same package-name-correct exports.
