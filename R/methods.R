#' @export
print.precfact <- function(x, ...) {
  ne <- sum(x$graph[upper.tri(x$graph)])
  tot <- x$d * (x$d - 1) / 2
  cat("Bayesian precision factor model (Chandra, Mueller and Sarkar 2026)\n")
  cat(sprintf("  n = %d observations, d = %d variables, q = %d factors\n", x$n, x$d, x$q))
  nch <- if (is.null(x$settings$chains)) 1L else x$settings$chains
  cat(sprintf("  %d chain(s) of %d MCMC iterations, %d burn-in, thinning %d; %.1f seconds\n",
              nch, x$settings$nrun, x$settings$burn, x$settings$thin, x$time))
  cat(sprintf("  Graph: %d edges of %d (%.2f%%); posterior FDR %.3f at epsilon = %.3f\n",
              ne, tot, 100 * ne / tot, x$pFDR, x$epsilon))
  if (!is.null(x$chain_agreement)) {
    a <- x$chain_agreement
    cat(sprintf("  Single chains select %s edges; each holds %.0f%%-%.0f%% of the pooled graph's edges\n",
                paste(a$n_edges, collapse = "/"), 100 * min(a$share_of_pooled_edges),
                100 * max(a$share_of_pooled_edges)))
  }
  invisible(x)
}

#' @export
summary.precfact <- function(object, top = 10, ...) {
  e <- precfact_edges(object)
  deg <- rowSums(object$graph)
  structure(list(n = object$n, d = object$d, q = object$q, n_edges = nrow(e),
                 density = nrow(e) / (object$d * (object$d - 1) / 2),
                 pFDR = object$pFDR, epsilon = object$epsilon,
                 positive = sum(e$partial_correlation > 0), negative = sum(e$partial_correlation < 0),
                 isolated = sum(deg == 0), hubs = utils::head(top_hubs(object), top),
                 strongest = utils::head(e, top)),
            class = "summary.precfact")
}

#' @export
print.summary.precfact <- function(x, ...) {
  cat(sprintf("Precision factor graph: %d variables, %d edges (density %.4f), q = %d\n",
              x$d, x$n_edges, x$density, x$q))
  cat(sprintf("  posterior FDR %.3f at epsilon = %.3f; %d positive and %d negative partial correlations; %d isolated variables\n",
              x$pFDR, x$epsilon, x$positive, x$negative, x$isolated))
  if (nrow(x$hubs)) {
    cat("\nMost connected variables:\n")
    print(x$hubs, row.names = FALSE, digits = 3)
    cat("\nStrongest edges:\n")
    print(x$strongest, row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' Edges of the selected graph
#'
#' Lists the edges of the graph selected by [precfact()], ordered by the absolute posterior
#' mean partial correlation.
#'
#' @param fit A `precfact` object.
#'
#' @return A data frame with `from`, `to`, `partial_correlation` (posterior mean), `sd`
#'   (posterior standard deviation) and `edge_prob` (posterior probability that
#'   \eqn{|\rho| > \epsilon} at the selected \eqn{\epsilon}).
#' @examples
#' sim <- simulate_precfact(n = 100, d = 15, q = 3, seed = 2)
#' fit <- precfact(sim$y, nrun = 500, burn = 150, seed = 2, verbose = FALSE)
#' head(precfact_edges(fit))
#' @export
precfact_edges <- function(fit) {
  .check_fit(fit)
  idx <- which(fit$graph & upper.tri(fit$graph), arr.ind = TRUE)
  e <- data.frame(from = fit$nodes[idx[, 1]], to = fit$nodes[idx[, 2]],
                  partial_correlation = fit$pcor[idx], sd = fit$pcor_sd[idx],
                  edge_prob = fit$edge_prob[idx], stringsAsFactors = FALSE)
  e[order(-abs(e$partial_correlation)), , drop = FALSE]
}

#' Most connected variables
#'
#' Ranks the variables of a fitted graph by their number of edges (degree), with ties
#' broken by strength, the sum of absolute partial correlations over their edges.
#'
#' @param fit A `precfact` object.
#' @param n Number of variables to return; `NULL` returns all with at least one edge.
#'
#' @return A data frame with `node`, `degree` and `strength`.
#' @examples
#' sim <- simulate_precfact(n = 100, d = 15, q = 3, seed = 2)
#' fit <- precfact(sim$y, nrun = 500, burn = 150, seed = 2, verbose = FALSE)
#' top_hubs(fit, 5)
#' @export
top_hubs <- function(fit, n = NULL) {
  .check_fit(fit)
  deg <- rowSums(fit$graph)
  str <- rowSums(abs(fit$pcor) * fit$graph)
  h <- data.frame(node = fit$nodes, degree = deg, strength = str, stringsAsFactors = FALSE,
                  row.names = NULL)
  h <- h[h$degree > 0, , drop = FALSE]
  h <- h[order(-h$degree, -h$strength), , drop = FALSE]
  rownames(h) <- NULL
  if (!is.null(n)) h <- utils::head(h, n)
  h
}

#' Selected partial correlation matrix
#'
#' Returns the posterior mean partial correlation matrix with the entries outside the
#' selected graph set to zero, optionally restricted to some variables.
#'
#' @param fit A `precfact` object.
#' @param nodes Optional names or indices of variables to keep.
#'
#' @return A symmetric matrix with unit diagonal.
#' @examples
#' sim <- simulate_precfact(n = 100, d = 12, q = 3, seed = 3)
#' fit <- precfact(sim$y, nrun = 500, burn = 150, seed = 3, verbose = FALSE)
#' round(graph_pcor(fit)[1:4, 1:4], 2)
#' @export
graph_pcor <- function(fit, nodes = NULL) {
  .check_fit(fit)
  m <- fit$pcor * fit$graph
  diag(m) <- 1
  if (!is.null(nodes)) m <- m[nodes, nodes, drop = FALSE]
  m
}

.check_fit <- function(fit) {
  if (!inherits(fit, "precfact")) stop("`fit` must be a precfact object.", call. = FALSE)
  invisible(TRUE)
}
