# Patch notes

## Added public input options

```r
ignore_branch_lengths = FALSE
ignore_node_sizes = FALSE
```

These are available in `all_indices()`, `node()`, `long_star()`,
`calculate_all_subtree_indices()`, and `process_tree_folder()`.
`read_convert()` accepts `ignore_branch_lengths`.

- Ignored branch lengths are replaced with unit edge lengths (`1`).
- Ignored node sizes cause `node_abundances` to be passed to C++ as `NULL`,
  which uses equal direct abundance on tips and zero direct abundance on
  internal nodes.
- Both options default to `FALSE`, preserving existing calls.

## Added package data

The patch includes two raw Newick files under `inst/extdata/`:

- `rphylo.nwk`: 500 tips, generated using `ape::rphylo(500, 1, 0)`;
- `sim-taxa.nwk`: 500 tips, generated using
  `TreeSimGM::sim.taxa(1, 500, waitsp = "rexp(1.2)")[[1]]`.

The exported `tree_index_example()` helper lists and locates these installed
files without relying on a machine-specific path.

## Added tests and examples

- `test-example-data.R` checks that both files are installed, readable, and
  have the expected 500-tip binary-tree structure.
- `test-basic-specification.R` checks the one-node boundary case,
  multifurcations, missing and negative branch lengths, and duplicate labels.
- The example-data tests run the full index pipeline only on ten-tip subsets,
  keeping routine checks fast.
- Executable roxygen examples were added to `read_convert()`, `all_indices()`,
  and `tree_index_example()`.

## Rcpp export repair

The included generated files use the package-name-correct
`_TreeShapeIndicesC_` symbol prefix. Running `Rcpp::compileAttributes()` from
the package root regenerates the same exports.
