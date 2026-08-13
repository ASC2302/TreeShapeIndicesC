test_that("a one-node tree is rejected", {
  one_node_tree <- structure(
    list(
      edge = matrix(integer(), nrow = 0L, ncol = 2L),
      tip.label = "A",
      Nnode = 0L
    ),
    class = "phylo"
  )

  expect_error(
    all_indices(one_node_tree),
    "no edges|at least two tips"
  )
})


test_that("a rooted multifurcating tree is accepted", {
  tree <- ape::read.tree(text = "(A:1,B:1,C:1,D:1);")

  expect_no_error(
    result <- all_indices(tree)
  )

  expect_named(
    result,
    c(
      "D0N", "D1N", "J1N",
      "D0S", "D1S", "J1S",
      "D0L", "D1L", "J1L"
    )
  )
})


test_that("missing branch lengths are replaced with unit lengths", {
  tree <- ape::read.tree(text = "((A,B),C);")
  converted <- read_convert(tree)

  expect_equal(
    converted$edge.length,
    rep(1, nrow(converted$edge))
  )
})


test_that("negative branch lengths are rejected unless ignored", {
  tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
  tree$edge.length[1] <- -1

  expect_error(
    all_indices(tree),
    "negative branch length"
  )

  expect_no_error(
    all_indices(tree, ignore_branch_lengths = TRUE)
  )
})


test_that("duplicate tip labels are rejected", {
  tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
  tree$tip.label[2] <- tree$tip.label[1]

  expect_error(
    all_indices(tree),
    "duplicate tip labels"
  )
})
