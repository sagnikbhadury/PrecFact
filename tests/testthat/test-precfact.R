test_that("fast inverse and determinant match base R", {
  set.seed(1)
  L <- matrix(rnorm(40), 10, 4)
  ps <- runif(10, 0.5, 2)
  Om <- tcrossprod(L) + diag(ps)
  expect_equal(fast_fact_inv(L, ps), solve(Om), tolerance = 1e-8)
  expect_equal(fast_fact_det(L, ps), as.numeric(determinant(Om)$modulus), tolerance = 1e-8)
})

test_that("pca_prec returns d x q with q >= 2", {
  set.seed(2)
  y <- matrix(rnorm(50 * 12), 50, 12)
  l <- pca_prec(y, nv = 8)
  expect_equal(nrow(l), 12)
  expect_gte(ncol(l), 2)
  expect_equal(ncol(pca_prec(y, nv = 8, q = 3)), 3)
  expect_equal(nrow(pca_prec(y, sp = FALSE)), 12)
})

test_that("sampler output has the documented shape and is reproducible", {
  y <- scale(simulate_precfact(n = 40, d = 8, q = 2, seed = 3)$y, scale = FALSE)
  li <- pca_prec(y, nv = 4)
  set.seed(10); r1 <- prec_est(Y = y, lambda_init = li, nrun = 60, thin = 3, verbose = FALSE)
  set.seed(10); r2 <- prec_est(Y = y, lambda_init = li, nrun = 60, thin = 3, verbose = FALSE)
  expect_equal(dim(r1$Lambda), c(8, ncol(li), 20))
  expect_equal(dim(r1$residuals), c(20, 8))
  expect_true(all(r1$residuals > 0))
  expect_identical(r1, r2)
  r3 <- prec_est(Y = y, lambda_init = li, nrun = 30, thin = 3, alldiff = TRUE, verbose = FALSE)
  expect_equal(dim(r3$residuals), c(10, 8))
})

test_that("post-processing agrees with a direct computation", {
  y <- scale(simulate_precfact(n = 40, d = 7, q = 2, seed = 4)$y, scale = FALSE)
  set.seed(4)
  res <- prec_est(Y = y, lambda_init = pca_prec(y, nv = 4), nrun = 100, thin = 2, verbose = FALSE)
  grid <- c(0.05, 0.1, 0.2)
  out <- fact_prec_processing(res, burnout = 10, epsilon_grid = grid)
  keep <- 11:50
  oms <- lapply(keep, function(s) tcrossprod(res$Lambda[, , s]) + diag(res$residuals[s, ]))
  expect_equal(out$Omega, Reduce(`+`, oms) / length(oms), tolerance = 1e-10)
  pc <- lapply(oms, function(o) { r <- -cov2cor(o); diag(r) <- 1; r })
  expect_equal(out$pcor, Reduce(`+`, pc) / length(pc), tolerance = 1e-10)
  g <- which(grid == out$epsilon)
  pp <- Reduce(`+`, lapply(oms, function(o) abs(cov2cor(o)) > grid[g])) / length(oms)
  diag(pp) <- 1
  expect_equal(out$edge_prob, pp, tolerance = 1e-12)
  expect_true(isSymmetric(out$P))
})

test_that("precfact recovers a clear graph and methods work", {
  sim <- simulate_precfact(n = 300, d = 15, q = 2, pr = 0.2, seed = 5)
  fit <- precfact(sim$y, nrun = 1500, burn = 500, seed = 5, verbose = FALSE)
  expect_s3_class(fit, "precfact")
  expect_true(isSymmetric(fit$graph))
  tg <- sim$graph[upper.tri(sim$graph)] & abs(sim$pcor[upper.tri(sim$pcor)]) > 0.1
  eg <- fit$graph[upper.tri(fit$graph)]
  expect_gt(sum(eg & tg) / max(sum(tg), 1), 0.5)
  expect_lt(sum(eg & !sim$graph[upper.tri(sim$graph)]) / max(sum(eg), 1), 0.3)
  e <- precfact_edges(fit)
  expect_equal(nrow(e), sum(eg))
  expect_s3_class(summary(fit), "summary.precfact")
  expect_output(print(fit), "precision factor")
  expect_s3_class(plot_network(fit), "ggplot")
  expect_s3_class(heatplot(fit), "ggplot")
  expect_s3_class(plot_trace(fit), "ggplot")
})

test_that("input errors are caught", {
  expect_error(precfact(matrix(1:6, 3, 2) + 0.5), "at least three")
  y <- matrix(rnorm(30), 10, 3); y[1, 1] <- NA
  expect_error(precfact(y), "missing")
  expect_error(precfact(matrix(rnorm(30), 10, 3), nrun = 10, burn = 20), "burn")
})

test_that("select_graph reproduces the fit and the default bounds the FDR", {
  sim <- simulate_precfact(n = 150, d = 12, q = 2, pr = 0.2, seed = 6)
  fit <- precfact(sim$y, nrun = 600, burn = 200, seed = 6, keep_samples = TRUE, verbose = FALSE)
  expect_lt(fit$pFDR, 0.1)
  same <- select_graph(fit)
  expect_identical(same$graph, fit$graph)
  loose <- suppressWarnings(select_graph(fit, alpha = 0.1, beta = 0.6))
  expect_gte(sum(loose$graph), 0)
  expect_error(select_graph(precfact(sim$y, nrun = 100, burn = 20, verbose = FALSE)), "keep_samples")
})
