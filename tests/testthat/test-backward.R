test_that("top_two returns the two closest columns as site ids", {
  sim = rbind(c(0.1, 0.9, 0.5),
               c(0.7, 0.2, 0.3))
  top = top_two(sim, cols = c(10L, 20L, 30L))
  expect_equal(top$best_site, c(20L, 10L))
  expect_equal(top$best, c(0.9, 0.7))
  expect_equal(top$second_site, c(30L, 30L))
  expect_equal(top$second, c(0.5, 0.3))
})

test_that("top_two has no runner-up with one column, and nothing with none", {
  one = top_two(matrix(c(0.4, 0.6), ncol = 1), cols = 5L)
  expect_equal(one$second_site, c(0L, 0L))
  expect_equal(one$second, c(0, 0))

  none = top_two(matrix(0, nrow = 2, ncol = 0), cols = integer(0))
  expect_equal(none$best_site, c(0L, 0L))
})

test_that("add_by_site accumulates by index and skips site 0", {
  total = add_by_site(c(1, 1, 1), site = c(2L, 2L, 3L, 0L), value = c(1, 2, 4, 100))
  expect_equal(total, c(1, 4, 5))
})

# Slow but obvious: coverage summed over species, each a mean over the sites
# holding it of the similarity to the closest retained holder.
coverage = function(kernel, presence, retained) {
  sum(vapply(seq_len(ncol(presence)), function(s) {
    holders = which(presence[, s] > 0)
    keep = intersect(holders, retained)
    if (length(keep) == 0) return(0)
    mean(apply(kernel[holders, keep, drop = FALSE], 1, max))
  }, numeric(1)))
}

naive_order = function(kernel, presence, weights) {
  retained = seq_len(nrow(kernel))
  order = integer(0)
  while (length(retained) > 0) {
    full = coverage(kernel, presence, retained)
    cost = vapply(retained, function(p) full - coverage(kernel, presence, setdiff(retained, p)),
                   numeric(1))
    p = retained[which.min(cost * weights[retained])]
    order = c(order, p)
    retained = setdiff(retained, p)
  }
  order
}

test_that("incremental elimination matches recomputing coverage from scratch", {
  for (seed in 1:5) {
    prob = toy_problem(n_sites = 14, n_species = 5, seed = seed)
    kernel = similarity(whiten(prob$env, q = 0.9), bw_mult = 1)
    weights = runif(14, 0.5, 3)
    res = backward_order(kernel, prob$species, min_sites = 0, weights = weights)
    expect_equal(res$order, naive_order(kernel, prob$species, weights), info = paste("seed", seed))
  }
})

test_that("the species floor protects sites until nothing else can go", {
  # species "rare" is held by sites 1 and 2 only; with a floor of 2 neither may
  # be dropped while any other site can be
  n = 8
  presence = cbind(common = rep(1, n), rare = c(1, 1, rep(0, n - 2)))
  set.seed(4)
  kernel = similarity(whiten(matrix(rnorm(n * 3), n), q = 1), 1)

  res = backward_order(kernel, presence, min_sites = 2, weights = rep(1, n))

  expect_setequal(res$order[(n - 1):n], c(1, 2))
  expect_equal(res$floor_released_at, 2L)
})
