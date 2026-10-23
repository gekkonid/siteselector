test_that("every site is ranked once, with the documented columns", {
  prob = toy_problem()
  out = rank_sites(prob$species, prob$env, min_sites = 0)

  expect_named(out, c("site", "rank", "loss", "weight"))
  expect_equal(out$rank, 1:12)
  expect_setequal(out$site, rownames(prob$species))
  expect_true(all(out$loss >= 0))
})

test_that("near-duplicate sites are not both picked early", {
  # two environmental clumps, all species everywhere: the best two sites are
  # one from each clump
  set.seed(5)
  ids = sprintf("s%d", 1:6)
  centres = rbind(matrix(0, 3, 2), matrix(10, 3, 2))
  env = centres + matrix(rnorm(12, sd = 0.1), 6, dimnames = list(ids, NULL))
  species = matrix(1, 6, 2, dimnames = list(ids, NULL))

  out = rank_sites(species, env, min_sites = 0)

  expect_equal(sort(as.integer(substr(out$site[1:2], 2, 2)) > 3), c(FALSE, TRUE))
})

test_that("the only holder of a species outranks a redundant neighbour", {
  ids = c("a", "b", "far1", "far2", "far3")
  env = rbind(a = c(0, 0), b = c(0.01, 0), far1 = c(5, 5), far2 = c(-5, 5), far3 = c(5, -5))
  species = cbind(common = rep(1, 5), only_b = c(0, 1, 0, 0, 0))
  rownames(species) = ids

  out = rank_sites(species, env, min_sites = 0)

  expect_lt(out$rank[out$site == "b"], out$rank[out$site == "a"])
})

test_that("a higher weight keeps an otherwise identical site longer", {
  prob = toy_problem(seed = 6)
  prob$env["s02", ] = prob$env["s01", ]
  prob$species["s02", ] = prob$species["s01", ]

  plain = rank_sites(prob$species, prob$env, min_sites = 0)
  favoured = rank_sites(prob$species, prob$env, weights = c(s02 = 3), min_sites = 0)

  expect_lt(favoured$rank[favoured$site == "s02"], favoured$rank[favoured$site == "s01"])
  expect_equal(favoured$weight[favoured$site == "s02"], 3)
  expect_true(all(plain$weight == 1))
})

test_that("rows of env may be in a different order from species", {
  prob = toy_problem()
  shuffled = prob$env[sample(nrow(prob$env)), ]

  expect_equal(rank_sites(prob$species, shuffled, min_sites = 0),
               rank_sites(prob$species, prob$env, min_sites = 0))
})

test_that("species held by no site are ignored", {
  prob = toy_problem()
  extra = cbind(prob$species, empty = 0)
  expect_equal(rank_sites(extra, prob$env, min_sites = 0),
               rank_sites(prob$species, prob$env, min_sites = 0))
})

test_that("a floor that cannot be kept gives a warning", {
  prob = toy_problem()
  expect_warning(rank_sites(prob$species, prob$env, min_sites = 100), "min_sites")
})

test_that("bad input is rejected", {
  prob = toy_problem()
  species = prob$species

  expect_error(rank_sites(species, prob$env[-1, ]), "same sites")
  bad_env = prob$env
  bad_env[1, 1] = NA
  expect_error(rank_sites(species, bad_env), "missing")
  expect_error(rank_sites(species, prob$env, weights = c(s01 = 0)), "positive")
  expect_error(rank_sites(species, prob$env, weights = c(nope = 2)), "unknown")
  expect_error(rank_sites(species, prob$env, weights = 2), "named")
  expect_error(rank_sites(species * 0, prob$env), "no species")
  unnamed = species
  rownames(unnamed) = NULL
  expect_error(rank_sites(unnamed, prob$env), "row names")
  expect_error(rank_sites(species, prob$env, q = 0), "`q`")
})
