# Bundled example trees

This directory contains the raw Newick example data installed with
`TreeShapeIndicesC`.

- `rphylo.nwk` contains a rooted, fully bifurcating tree with 500 tips,
  generated with `ape::rphylo(500, 1, 0)`.
- `sim-taxa.nwk` contains a rooted, fully bifurcating tree with 500 tips,
  generated with `TreeSimGM::sim.taxa(1, 500, waitsp = "rexp(1.2)")[[1]]`.

Both files have unique tip labels and one finite, positive branch length for
each edge. They are deliberately stored as raw files under `inst/extdata/` so
that the package's file-reading interface can be demonstrated and tested.

Use `tree_index_example()` to locate the installed copies. Do not construct a
path to `inst/extdata/` directly, because `inst/` is removed when an R package
is installed.
