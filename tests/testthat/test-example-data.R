test_that("bundled example tree files are installed", {
  expect_setequal(
    tree_index_example(),
    c("rphylo.nwk", "sim-taxa.nwk")
  )

  paths <- vapply(
    tree_index_example(),
    tree_index_example,
    character(1)
  )

  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size > 0))
})


test_that("tree_index_example rejects unknown filenames", {
  expect_error(
    tree_index_example("unknown-tree.nwk"),
    "must be NULL or one of"
  )

  expect_error(
    tree_index_example(1),
    "must be NULL or one of"
  )
})


test_that("rphylo example data have the expected structure", {
  tree <- read_convert(tree_index_example("rphylo.nwk"))

  expect_s3_class(tree, "phylo")
  expect_equal(ape::Ntip(tree), 500L)
  expect_equal(tree$Nnode, 499L)
  expect_equal(nrow(tree$edge), 998L)
  expect_false(anyDuplicated(tree$tip.label) > 0L)
  expect_length(tree$edge.length, nrow(tree$edge))
  expect_true(all(is.finite(tree$edge.length)))
  expect_true(all(tree$edge.length > 0))
})


test_that("sim-taxa example data have the expected structure", {
  tree <- read_convert(tree_index_example("sim-taxa.nwk"))

  expect_s3_class(tree, "phylo")
  expect_equal(ape::Ntip(tree), 500L)
  expect_equal(tree$Nnode, 499L)
  expect_equal(nrow(tree$edge), 998L)
  expect_false(anyDuplicated(tree$tip.label) > 0L)
  expect_length(tree$edge.length, nrow(tree$edge))
  expect_true(all(is.finite(tree$edge.length)))
  expect_true(all(tree$edge.length > 0))
})


test_that("bundled trees can be used in a quick index calculation", {
  expected_names <- c(
    "D0N", "D1N", "J1N",
    "D0S", "D1S", "J1S",
    "D0L", "D1L", "J1L"
  )

  for (filename in tree_index_example()) {
    tree <- read_convert(tree_index_example(filename))
    small_tree <- ape::keep.tip(tree, tree$tip.label[seq_len(10)])

    result <- all_indices(
      small_tree,
      ignore_branch_lengths = TRUE
    )

    expect_named(result, expected_names)

    values <- vapply(result[expected_names], as.numeric, numeric(1))
    expect_true(all(is.finite(values)))
    expect_true(all(values[c("D0N", "D1N", "D0S", "D1S", "D0L", "D1L")] > 0))
    expect_true(all(values[c("J1N", "J1S", "J1L")] >= -1e-8))
    expect_true(all(values[c("J1N", "J1S", "J1L")] <= 1 + 1e-8))
  }
})
