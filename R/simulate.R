#' Simulate data with a sparse low-rank-plus-diagonal precision matrix
#'
#' Generates \eqn{y_i \sim N_d(0, \Omega^{-1})} with \eqn{\Omega = \Lambda\Lambda^T + \psi I}
#' and a sparse \eqn{\Lambda}: each of its `q` columns has a run of
#' `floor(d * pr)` consecutive non-zero entries drawn from \eqn{U(-val, val)} at a random
#' position. This is the design of the example in the original PrecFact vignette
#' (`simulate_lambda_sparse()`).
#'
#' @param n Sample size.
#' @param d Dimension.
#' @param q Number of columns of \eqn{\Lambda}.
#' @param pr Fraction of non-zero entries in each column.
#' @param val Range of the non-zero loadings.
#' @param psi Diagonal term.
#' @param seed Optional seed.
#'
#' @return A list with `y` (\eqn{n \times d}), `Omega`, `Lambda`, `pcor` (true partial
#'   correlations) and `graph` (true adjacency, \eqn{|\rho_{ij}| > 0}).
#' @examples
#' sim <- simulate_precfact(n = 50, d = 20, q = 3, seed = 1)
#' dim(sim$y)
#' @export
simulate_precfact <- function(n = 200, d = 50, q = 5, pr = 0.1, val = 2.5, psi = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- max(floor(d * pr), 1)
  Lambda <- vapply(seq_len(q), function(h) {
    x <- numeric(d)
    j <- sample.int(d - k, 1)
    x[seq_len(k) + j] <- stats::runif(k, -val, val)
    x
  }, numeric(d))
  Lambda <- matrix(Lambda, d, q)
  Omega <- tcrossprod(Lambda) + diag(psi, d)
  R <- chol(Omega)
  z <- matrix(stats::rnorm(n * d), n, d)
  y <- t(backsolve(R, t(z)))
  nm <- paste0("V", seq_len(d))
  colnames(y) <- nm
  dimnames(Omega) <- list(nm, nm)
  isd <- 1 / sqrt(diag(Omega))
  pc <- -Omega * outer(isd, isd)
  diag(pc) <- 1
  graph <- abs(pc) > 1e-10
  diag(graph) <- FALSE
  list(y = y, Omega = Omega, Lambda = Lambda, pcor = pc, graph = graph)
}
