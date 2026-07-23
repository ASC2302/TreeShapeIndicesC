test_that("ignored branch lengths are replaced with unit edges", {
  tree <- ape::read.tree(text = "((A:2,B:4):3,C:5);")

  converted <- read_convert(tree, ignore_branch_lengths = TRUE)

  expect_equal(converted$edge.length, rep(1, nrow(converted$edge)))
})


test_that("invalid lengths can be deliberately ignored", {
  tree <- ape::read.tree(text = "((A:2,B:4):3,C:5);")
  tree$edge.length[1] <- -2

  expect_error(read_convert(tree), "negative branch length")
  expect_silent(read_convert(tree, ignore_branch_lengths = TRUE))
})


test_that("topology-only output equals an explicitly unit-length tree", {
  tree <- ape::read.tree(text = "((A:2,B:4):3,C:5);")
  unit_tree <- tree
  unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))

  ignored <- all_indices(tree, ignore_branch_lengths = TRUE)
  explicit <- all_indices(unit_tree)

  expect_equal(ignored, explicit, tolerance = 1e-10)
})


test_that("ignored node sizes equal the default abundance calculation", {
  tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
  node_sizes <- data.frame(
    names = c("A", "B", "C"),
    values = c(0.8, 0.1, 0.1),
    stringsAsFactors = FALSE
  )

  ignored <- all_indices(
    tree,
    node_abundances = node_sizes,
    ignore_node_sizes = TRUE
  )
  default <- all_indices(tree)

  expect_equal(ignored, default, tolerance = 1e-10)
})


test_that("logical switches reject invalid values", {
  tree <- ape::read.tree(text = "(A:1,B:1);")

  expect_error(
    all_indices(tree, ignore_branch_lengths = NA),
    "must be either TRUE or FALSE"
  )
  expect_error(
    all_indices(tree, ignore_node_sizes = 1),
    "must be either TRUE or FALSE"
  )
})
