# Sampler interface. Method: Chandra, Mueller and Sarkar (2026), Bayesian Analysis,
# doi:10.1214/24-BA1461. Original code: Noirrit Kiran Chandra.

#' PCA-based initial value and latent dimension
#'
#' Chooses the number of columns \eqn{q} of \eqn{\Lambda} and an initial value for it from
#' a principal component analysis of the data, as in Section 2.2 of Chandra et al. (2026):
#' the inverse singular values, in decreasing order, are accumulated until they reach the
#' fraction `pr` of their total, and \eqn{q} is the number needed (at least 2). The initial
#' loadings are the corresponding right singular vectors scaled by
#' \eqn{\sqrt{n}/s_h}, where \eqn{s_h} are the singular values.
#'
#' @param Y An \eqn{n \times d} numeric data matrix, one observation per row.
#' @param nv Number of singular values to compute with `irlba` when `sp = TRUE`; an upper
#'   bound on \eqn{q}.
#' @param pr Fraction of the sum of inverse singular values to retain. The article uses
#'   0.95.
#' @param sp Use the truncated singular value decomposition of `irlba::irlba()`
#'   (Baglama and Reichel 2005)? If `FALSE`, the full `svd()` is used.
#' @param q Optional fixed number of columns, overriding `pr`.
#'
#' @return A \eqn{d \times q} matrix, the initial value of \eqn{\Lambda}.
#' @references
#' Chandra NK, Mueller P, Sarkar A (2026). Bayesian Scalable Precision Factor Analysis for
#' Gaussian Graphical Models. Bayesian Analysis, 21(1), 75-103. \doi{10.1214/24-BA1461}
#'
#' Baglama J, Reichel L (2005). Augmented implicitly restarted Lanczos bidiagonalization
#' methods. SIAM Journal on Scientific Computing, 27, 19-42. \doi{10.1137/04060593X}
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(100 * 20), 100, 20)
#' dim(pca_prec(y, nv = 10))
#' @export
pca_prec <- function(Y, nv = 25, pr = 0.95, sp = TRUE, q = NULL) {
  Y <- .check_data(Y)
  if (sp) {
    nv <- min(nv, nrow(Y) - 1, ncol(Y) - 1)
    if (nv > min(dim(Y)) / 2) {
      # irlba is unreliable for a large share of the singular values; the exact SVD,
      # truncated to the same nv leading values, gives the same quantity.
      pca <- svd(Y, nu = 0, nv = nv)
      pca$d <- pca$d[seq_len(nv)]
    } else {
      pca <- irlba::irlba(Y, nv = nv)
    }
  } else {
    pca <- svd(Y, nu = 0)
  }
  keep <- pca$d > max(pca$d) * 1e-10
  dvals <- pca$d[keep]
  V <- pca$v[, keep, drop = FALSE]
  nq <- length(dvals)
  if (is.null(q)) {
    inv <- 1 / dvals[nq:1]
    cum <- cumsum(inv) / sum(inv)
    q <- max(min(which(cum > pr)), 2)
  }
  q <- min(as.integer(q), nq)
  cols <- (1 + nq - q):nq
  V[, cols, drop = FALSE] %*% diag(1 / dvals[cols], nrow = q) * sqrt(nrow(Y))
}

#' Gibbs sampler for the precision factor model
#'
#' Draws from the posterior of \eqn{\Lambda} and \eqn{\Delta} in
#' \eqn{y_i \sim N_d(0, \Omega^{-1})}, \eqn{\Omega = \Lambda\Lambda^T + \Delta}, with the
#' Gibbs sampler of Section 2.3 of Chandra et al. (2026). Most users should call
#' [precfact()], which also chooses the initial value and selects the graph.
#'
#' The prior is \eqn{vec(\Lambda) \sim DL(a, b)}, the two-parameter Dirichlet-Laplace prior
#' (Bhattacharya et al. 2015), and \eqn{\delta_j^2 \sim G}, \eqn{G \sim DP(M, Ga(a_s, b_s))},
#' \eqn{M \sim Ga(a_{dir}, b_{dir})}. With `alldiff = TRUE` the \eqn{\delta_j^2} are instead
#' independent \eqn{Ga(a_s, b_s)}.
#'
#' @param as,bs Shape and rate of the gamma base measure of the Dirichlet process on the
#'   diagonal of \eqn{\Delta}.
#' @param a,b Parameters of the Dirichlet-Laplace prior on \eqn{\Lambda}.
#' @param nrun Number of MCMC iterations.
#' @param thin Thinning interval; `floor(nrun / thin)` draws are saved.
#' @param lambda_init Initial \eqn{d \times q} loading matrix, for example from
#'   [pca_prec()]. Its number of columns fixes \eqn{q}.
#' @param Y \eqn{n \times d} data matrix with centred columns.
#' @param alldiff If `TRUE`, the diagonal entries of \eqn{\Delta} are independent instead
#'   of clustered by a Dirichlet process.
#' @param a_dir,b_dir Shape and rate of the gamma prior on the Dirichlet process
#'   concentration.
#' @param verbose Print progress?
#'
#' @return A list with `Lambda`, a \eqn{d \times q \times S} array of draws of
#'   \eqn{\Lambda}; `residuals`, an \eqn{S \times d} matrix of draws of the diagonal of
#'   \eqn{\Delta}; and `n_clusters`, the number of distinct diagonal values in each draw.
#' @references
#' Chandra NK, Mueller P, Sarkar A (2026). Bayesian Scalable Precision Factor Analysis for
#' Gaussian Graphical Models. Bayesian Analysis, 21(1), 75-103. \doi{10.1214/24-BA1461}
#'
#' Bhattacharya A, Pati D, Pillai NS, Dunson DB (2015). Dirichlet-Laplace priors for optimal
#' shrinkage. Journal of the American Statistical Association, 110, 1479-1490.
#' \doi{10.1080/01621459.2014.960967}
#' @examples
#' set.seed(1)
#' y <- scale(matrix(rnorm(60 * 10), 60, 10), scale = FALSE)
#' res <- prec_est(Y = y, lambda_init = pca_prec(y, nv = 5), nrun = 100, thin = 5,
#'                 verbose = FALSE)
#' dim(res$Lambda)
#' @export
prec_est <- function(as = 0.1, bs = 0.1, a = 0.5, b = 0.5, nrun = 5500, thin = 5,
                     lambda_init, Y, alldiff = FALSE, a_dir = 0.1, b_dir = 0.1,
                     verbose = TRUE) {
  Y <- .check_data(Y)
  if (!is.matrix(lambda_init) || nrow(lambda_init) != ncol(Y))
    stop("`lambda_init` must be a matrix with one row per column of `Y`.", call. = FALSE)
  .check_positive(as = as, bs = bs, a = a, b = b, a_dir = a_dir, b_dir = b_dir)
  nrun <- as.integer(nrun)
  thin <- as.integer(thin)
  if (is.na(thin) || thin < 1L) stop("`thin` must be a positive integer.", call. = FALSE)
  if (is.na(nrun) || nrun < thin) stop("`nrun` must be at least `thin`.", call. = FALSE)
  prec_est_cpp(as, bs, a, b, nrun, thin, lambda_init, Y, isTRUE(alldiff), a_dir, b_dir,
               isTRUE(verbose))
}

#' Fit the Bayesian precision factor model and select the graph
#'
#' Fits the Gaussian precision factor model of Chandra, Mueller and Sarkar (2026) and
#' selects the conditional independence graph by posterior false discovery rate control.
#' The steps are: centre (and by default scale) the columns of `Y`; choose \eqn{q} and an
#' initial \eqn{\Lambda} with [pca_prec()]; run the Gibbs sampler [prec_est()]; discard
#' the burn-in; and apply the edge selection rule of Section 2.4 of the article (see
#' [fact_prec_processing()]).
#'
#' The defaults are the settings used for every analysis in the article: \eqn{a = b = 0.5},
#' \eqn{a_s = b_s = a_{dir} = b_{dir} = 0.1}, 5,500 iterations with the first 1,250
#' discarded and a thinning interval of 5, \eqn{q} chosen to cover 95% of the inverse
#' singular values, and the posterior FDR controlled at 0.10 (\eqn{\beta = 0.9}).
#'
#' @param Y An \eqn{n \times d} numeric matrix or data frame, one observation per row.
#'   Missing values are not allowed.
#' @param q Number of columns of \eqn{\Lambda}. By default chosen from the data by
#'   [pca_prec()].
#' @param nv,pr,sparse_pca Passed to [pca_prec()] as `nv`, `pr` and `sp`.
#' @param as,bs,a,b,a_dir,b_dir,alldiff Prior settings; see [prec_est()].
#' @param nrun Number of MCMC iterations.
#' @param burn Number of initial iterations to discard. Must be smaller than `nrun`.
#' @param thin Thinning interval.
#' @param alpha Target posterior false discovery rate.
#' @param beta Threshold on the posterior probability of an edge; by default
#'   \eqn{1 - \alpha}, which bounds the posterior FDR by \eqn{\alpha} (Section 2.4 of the
#'   article). See [fact_prec_processing()].
#' @param epsilon_grid Grid of values for the interval null \eqn{|\rho_{ij}| \le \epsilon}.
#' @param center,scale Centre and scale the columns of `Y` before fitting? The model has
#'   mean zero, so centring should only be turned off for data that are already centred.
#'   With `scale = TRUE` the precision matrix refers to the standardised variables; the
#'   partial correlations and the graph do not depend on the scale.
#' @param chains Number of independent chains, run one after another from the same
#'   initial value. Their post-burn-in draws are pooled for every posterior summary and for
#'   the graph, so an edge is selected only if it is supported across chains. With more
#'   than one chain, `chain_agreement` compares the graph each chain gives on its own.
#'   Several chains are recommended for real data, where single chains can settle on
#'   different graphs (see the README).
#' @param seed Optional seed passed to [set.seed()] before sampling.
#' @param keep_samples Keep the pooled post-burn-in draws of \eqn{\Lambda} and
#'   \eqn{\Delta} in the result? They take \eqn{(d q + d) S} numbers.
#' @param verbose Print progress?
#'
#' @return An object of class `precfact`, a list with
#' \describe{
#'   \item{`Omega`}{Posterior mean of the precision matrix.}
#'   \item{`pcor`, `pcor_sd`}{Posterior mean and standard deviation of the partial
#'     correlation matrix \eqn{\rho_{ij} = -\omega_{ij}/(\omega_{ii}\omega_{jj})^{1/2}}.}
#'   \item{`graph`}{Logical adjacency matrix of the selected graph.}
#'   \item{`edge_prob`}{Posterior probabilities \eqn{\Pi(|\rho_{ij}| > \epsilon \mid y)} at
#'     the selected \eqn{\epsilon}.}
#'   \item{`pFDR`, `epsilon`}{Posterior FDR of the selected graph and the selected
#'     \eqn{\epsilon}.}
#'   \item{`fdr_table`}{Number of edges and posterior FDR for every value of the grid.}
#'   \item{`logdet`, `chain`}{Log-determinant of \eqn{\Omega} in each kept draw, and the
#'     chain it came from, for checking convergence.}
#'   \item{`chain_agreement`}{With `chains > 1`: a data frame with one row per chain giving
#'     the number of edges that chain selects alone, the share of the pooled graph's edges
#'     it contains, and the correlation of its posterior mean partial correlations with
#'     those of the other chains (average over pairs); `NULL` for one chain.}
#'   \item{`n_clusters`}{Number of distinct diagonal values of \eqn{\Delta} in each saved
#'     draw.}
#'   \item{`q`, `n`, `d`, `nodes`, `settings`, `time`}{Dimensions, variable names,
#'     settings and run time in seconds.}
#'   \item{`samples`}{The posterior draws, if `keep_samples = TRUE`.}
#' }
#' @references
#' Chandra NK, Mueller P, Sarkar A (2026). Bayesian Scalable Precision Factor Analysis for
#' Gaussian Graphical Models. Bayesian Analysis, 21(1), 75-103. \doi{10.1214/24-BA1461}
#' @seealso [plot_network()], [heatplot()], [precfact_edges()], [top_hubs()]
#' @examples
#' sim <- simulate_precfact(n = 100, d = 20, q = 3, seed = 1)
#' fit <- precfact(sim$y, nrun = 600, burn = 200, seed = 1, verbose = FALSE)
#' fit
#' head(precfact_edges(fit))
#' @export
precfact <- function(Y, q = NULL, nv = 25, pr = 0.95, sparse_pca = TRUE,
                     as = 0.1, bs = 0.1, a = 0.5, b = 0.5, a_dir = 0.1, b_dir = 0.1,
                     alldiff = FALSE, nrun = 5500, burn = 1250, thin = 5,
                     alpha = 0.1, beta = 1 - alpha, epsilon_grid = seq(0.05, 0.3, length.out = 15),
                     chains = 1, center = TRUE, scale = TRUE, seed = NULL,
                     keep_samples = FALSE, verbose = TRUE) {
  cl <- match.call()
  Y <- .check_data(Y)
  if (burn >= nrun) stop("`burn` must be smaller than `nrun`.", call. = FALSE)
  if (ncol(Y) < 3) stop("`Y` must have at least three columns.", call. = FALSE)
  nodes <- colnames(Y)
  if (is.null(nodes)) nodes <- paste0("V", seq_len(ncol(Y)))
  Ys <- base::scale(Y, center = center, scale = scale)
  if (any(!is.finite(Ys))) stop("A column of `Y` is constant; remove it before fitting.", call. = FALSE)
  attributes(Ys) <- list(dim = dim(Y))
  chains <- as.integer(chains)
  if (is.na(chains) || chains < 1L) stop("`chains` must be a positive integer.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  lam_init <- pca_prec(Ys, nv = nv, pr = pr, sp = sparse_pca, q = q)
  if (verbose) message(sprintf("PrecFact: n = %d, d = %d, q = %d, %d chain(s)", nrow(Ys), ncol(Ys),
                               ncol(lam_init), chains))
  burnout <- floor(burn / thin)
  runs <- lapply(seq_len(chains), function(ch) {
    if (verbose && chains > 1) message("PrecFact: chain ", ch)
    r <- prec_est(as = as, bs = bs, a = a, b = b, nrun = nrun, thin = thin,
                  lambda_init = lam_init, Y = Ys, alldiff = alldiff, a_dir = a_dir,
                  b_dir = b_dir, verbose = verbose)
    keep <- seq.int(burnout + 1L, dim(r$Lambda)[3])
    list(Lambda = r$Lambda[, , keep, drop = FALSE], residuals = r$residuals[keep, , drop = FALSE],
         n_clusters = r$n_clusters[keep])
  })
  # pooled post-burn-in draws of all chains
  res <- list(Lambda = array(unlist(lapply(runs, `[[`, "Lambda")),
                             c(ncol(Ys), ncol(lam_init), sum(vapply(runs, function(r) dim(r$Lambda)[3], 1)))),
              residuals = do.call(rbind, lapply(runs, `[[`, "residuals")),
              n_clusters = unlist(lapply(runs, `[[`, "n_clusters")))
  out <- fact_prec_processing(res, burnout = 0, alpha = alpha, beta = beta,
                              epsilon_grid = epsilon_grid)
  chain_id <- rep(seq_len(chains), vapply(runs, function(r) dim(r$Lambda)[3], 1))
  agreement <- NULL
  if (chains > 1) {
    per <- lapply(runs, function(r) suppressWarnings(
      fact_prec_processing(r, burnout = 0, alpha = alpha, beta = beta, epsilon_grid = epsilon_grid)))
    agreement <- .chain_agreement(per, out)
  }
  dn <- list(nodes, nodes)
  fit <- list(Omega = out$Omega, pcor = out$pcor, pcor_sd = out$pcor_sd,
              graph = out$P == 1, edge_prob = out$edge_prob,
              pFDR = out$pFDR, epsilon = out$epsilon, fdr_table = out$fdr_table,
              logdet = out$logdet, chain = chain_id, n_clusters = res$n_clusters,
              chain_agreement = agreement,
              q = ncol(lam_init), n = nrow(Ys), d = ncol(Ys), nodes = nodes,
              settings = list(as = as, bs = bs, a = a, b = b, a_dir = a_dir, b_dir = b_dir,
                              alldiff = alldiff, nrun = nrun, burn = burn, thin = thin,
                              chains = chains, alpha = alpha, beta = beta,
                              epsilon_grid = epsilon_grid, center = center, scale = scale,
                              seed = seed),
              time = proc.time()[["elapsed"]] - t0, call = cl)
  for (nm in c("Omega", "pcor", "pcor_sd", "graph", "edge_prob")) dimnames(fit[[nm]]) <- dn
  if (keep_samples) fit$samples <- res
  class(fit) <- "precfact"
  fit
}

.check_data <- function(Y) {
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y)) stop("`Y` must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(Y)) stop("`Y` contains missing values.", call. = FALSE)
  storage.mode(Y) <- "double"
  Y
}

.check_positive <- function(...) {
  x <- list(...)
  bad <- names(x)[!vapply(x, function(v) is.numeric(v) && length(v) == 1 && is.finite(v) && v > 0, logical(1))]
  if (length(bad)) stop("These settings must be positive numbers: ", paste(bad, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.chain_agreement <- function(per, pooled) {
  up <- upper.tri(pooled$P)
  pg <- pooled$P[up] == 1
  pc <- lapply(per, function(o) o$pcor[up])
  k <- length(per)
  data.frame(
    chain = seq_len(k),
    n_edges = vapply(per, function(o) sum(o$P[up]), numeric(1)),
    share_of_pooled_edges = vapply(per, function(o) if (any(pg)) mean(o$P[up][pg] == 1) else NA_real_, numeric(1)),
    pcor_cor_with_others = vapply(seq_len(k), function(i)
      mean(vapply(setdiff(seq_len(k), i), function(j) stats::cor(pc[[i]], pc[[j]]), numeric(1))), numeric(1)))
}
