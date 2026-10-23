test_that("whitened distance is Mahalanobis distance", {
  set.seed(2)
  # strongly correlated variables, which plain scaling would double count
  x = rnorm(40)
  env = cbind(a = x, b = x + rnorm(40, sd = 0.3), c = rnorm(40))

  z = whiten(env, q = 1)

  expected = sqrt(stats::mahalanobis(env[2, ], env[1, ], stats::cov(env)))
  expect_equal(sqrt(sum((z[2, ] - z[1, ])^2)), expected)
})

test_that("whitening keeps the axes needed for the requested variance", {
  set.seed(3)
  env = cbind(rnorm(50), rnorm(50), rnorm(50), rnorm(50))
  expect_equal(ncol(whiten(env, q = 1)), 4)
  expect_lt(ncol(whiten(env, q = 0.5)), 4)
})

test_that("constant variables are ignored rather than breaking scaling", {
  env = cbind(varying = c(1, 2, 3, 5, 8), constant = 4)
  expect_equal(ncol(whiten(env, q = 1)), 1)
})

test_that("similarity is a symmetric kernel with ones on the diagonal", {
  z = whiten(toy_problem()$env, q = 1)
  k = similarity(z, bw_mult = 1)
  expect_equal(unname(diag(k)), rep(1, nrow(z)))
  expect_equal(k, t(k))
  expect_true(all(k > 0 & k <= 1))
})

test_that("a wider bandwidth makes distinct sites more similar", {
  z = whiten(toy_problem()$env, q = 1)
  expect_gt(similarity(z, 3)[1, 2], similarity(z, 1)[1, 2])
})
