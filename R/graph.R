# Edge selection by posterior false discovery rate control, Section 2.4 of
# Chandra, Mueller and Sarkar (2026), doi:10.1214/24-BA1461.

#' Posterior summaries and graph selection from the MCMC draws
#'
#' Computes the posterior mean of \eqn{\Omega} and of the partial correlation matrix, and
#' selects the graph. For each \eqn{\epsilon} in `epsilon_grid` the edge \eqn{(i, j)} is
#' declared when \eqn{\Pi(|\rho_{ij}| > \epsilon \mid y) > \beta}, and the posterior false
#' discovery rate of these decisions is
#' \deqn{FDR(\epsilon) = \sum_{i<j} d_{ij}\{1 - \Pi(|\rho_{ij}| > \epsilon \mid y)\} /
#' \max(\sum_{i<j} d_{ij}, 1).}
#' Among the values of \eqn{\epsilon} with \eqn{FDR(\epsilon) < \alpha}, the one with the
#' largest FDR is selected, as in the original code. If none reaches \eqn{\alpha}, the one
#' with the smallest FDR is used with a warning.
#'
#' Because every declared edge has posterior probability above \eqn{\beta}, the posterior
#' FDR is at most \eqn{1 - \beta}. The article controls it at 0.10 with this bound, so the
#' default is \eqn{\beta = 1 - \alpha = 0.9}. (The original code used \eqn{\beta = 0.7}, for
#' which no \eqn{\epsilon} in the grid may reach the target; pass `beta = 0.7` to reproduce
#' it.)
#'
#' The draws are processed one at a time, so memory grows with \eqn{d^2} times the length
#' of the grid rather than with \eqn{d^2} times the number of draws.
#'
#' @param res Output of [prec_est()].
#' @param burnout Number of saved draws (after thinning) to discard.
#' @param alpha Target posterior false discovery rate.
#' @param beta Threshold on the posterior edge probability.
#' @param epsilon_grid Grid of values of \eqn{\epsilon}.
#'
#' @return A list with `Omega` (posterior mean), `P` (0/1 adjacency matrix of the selected
#'   graph), `pFDR`, and also `pcor`, `pcor_sd`, `edge_prob`, `epsilon`, `fdr_table` and
#'   `logdet` (see [precfact()]).
#' @references
#' Chandra NK, Mueller P, Sarkar A (2026). Bayesian Scalable Precision Factor Analysis for
#' Gaussian Graphical Models. Bayesian Analysis, 21(1), 75-103. \doi{10.1214/24-BA1461}
#'
#' Mueller P, Parmigiani G, Robert C, Rousseau J (2004). Optimal sample size for multiple
#' testing: the case of gene expression microarrays. Journal of the American Statistical
#' Association, 99, 990-1001.
#' @examples
#' set.seed(1)
#' y <- scale(matrix(rnorm(60 * 8), 60, 8), scale = FALSE)
#' res <- prec_est(Y = y, lambda_init = pca_prec(y, nv = 4), nrun = 200, thin = 5,
#'                 verbose = FALSE)
#' out <- fact_prec_processing(res, burnout = 10)
#' out$pFDR
#' @export
fact_prec_processing <- function(res, burnout = 100, alpha = 0.1, beta = 1 - alpha,
                                 epsilon_grid = seq(0.05, 0.3, length.out = 15)) {
  S <- dim(res$Lambda)[3]
  if (burnout >= S) stop("`burnout` must be smaller than the number of saved draws (", S, ").", call. = FALSE)
  if (!is.numeric(epsilon_grid) || any(epsilon_grid <= 0 | epsilon_grid >= 1))
    stop("`epsilon_grid` must lie in (0, 1).", call. = FALSE)
  keep <- seq.int(burnout + 1L, S) - 1L
  sm <- posterior_summaries_cpp(res$Lambda, res$residuals, keep, as.numeric(epsilon_grid))
  up <- upper.tri(sm$Omega)
  fdr <- vapply(seq_along(epsilon_grid), function(g) {
    pp <- sm$postprob[, , g][up]
    dec <- pp > beta
    c(sum(dec), sum(dec * (1 - pp)) / max(sum(dec), 1))
  }, numeric(2))
  fdr_table <- data.frame(epsilon = epsilon_grid, n_edges = fdr[1, ], pFDR = fdr[2, ])
  ok <- which(fdr_table$pFDR < alpha)
  if (length(ok)) {
    ind <- ok[which.max(fdr_table$pFDR[ok])]
  } else {
    warning("No value of epsilon reaches a posterior FDR below alpha = ", alpha,
            "; using the smallest FDR (", signif(min(fdr_table$pFDR), 3),
            "). Consider a larger beta.", call. = FALSE)
    ind <- which.min(fdr_table$pFDR)
  }
  edge_prob <- sm$postprob[, , ind]
  P <- (edge_prob > beta) * 1
  diag(P) <- 0
  diag(edge_prob) <- 1
  list(Omega = sm$Omega, P = P, pFDR = fdr_table$pFDR[ind], pcor = sm$pcor,
       pcor_sd = sm$pcor_sd, edge_prob = edge_prob, epsilon = epsilon_grid[ind],
       fdr_table = fdr_table, logdet = sm$logdet)
}

#' Re-select the graph of a fit with other settings
#'
#' Applies the edge selection of [fact_prec_processing()] again with a different target
#' FDR, threshold or grid, without re-running the sampler. The fit must have been made with
#' `keep_samples = TRUE`.
#'
#' @param fit A `precfact` object fitted with `keep_samples = TRUE`.
#' @inheritParams fact_prec_processing
#'
#' @return The `precfact` object with the graph, posterior summaries and settings updated.
#' @examples
#' sim <- simulate_precfact(n = 100, d = 15, q = 3, seed = 1)
#' fit <- precfact(sim$y, nrun = 500, burn = 150, seed = 1, keep_samples = TRUE,
#'                 verbose = FALSE)
#' fit2 <- select_graph(fit, alpha = 0.05)
#' c(sum(fit$graph), sum(fit2$graph)) / 2
#' @export
select_graph <- function(fit, alpha = 0.1, beta = 1 - alpha,
                         epsilon_grid = fit$settings$epsilon_grid) {
  .check_fit(fit)
  if (is.null(fit$samples))
    stop("The fit has no posterior draws; refit with `keep_samples = TRUE`.", call. = FALSE)
  out <- fact_prec_processing(fit$samples, burnout = 0, alpha = alpha,
                              beta = beta, epsilon_grid = epsilon_grid)
  dn <- list(fit$nodes, fit$nodes)
  fit$graph <- out$P == 1
  fit$edge_prob <- out$edge_prob
  dimnames(fit$graph) <- dimnames(fit$edge_prob) <- dn
  fit$pFDR <- out$pFDR
  fit$epsilon <- out$epsilon
  fit$fdr_table <- out$fdr_table
  fit$settings[c("alpha", "beta", "epsilon_grid")] <- list(alpha, beta, epsilon_grid)
  fit
}
