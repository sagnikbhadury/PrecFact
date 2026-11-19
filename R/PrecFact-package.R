#' PrecFact: Bayesian Precision Factor Analysis for Gaussian Graphical Models
#'
#' Estimates a high-dimensional Gaussian precision matrix and its conditional independence
#' graph with the precision factor model of Chandra, Mueller and Sarkar (2026). The
#' precision matrix is decomposed as \eqn{\Omega = \Lambda\Lambda^T + \Delta}, with a
#' \eqn{d \times q} loading matrix \eqn{\Lambda} under a Dirichlet-Laplace prior and a
#' diagonal \eqn{\Delta} whose entries follow a Dirichlet process. A latent variable
#' representation of this decomposition gives a Gibbs sampler that scales to thousands of
#' variables, and edges are selected by controlling the posterior false discovery rate.
#'
#' The main function is [precfact()]. [plot_network()] and [heatplot()] draw the estimated
#' graph and partial correlation matrix; [precfact_edges()] and [top_hubs()] list edges and
#' hub variables. The low-level sampler [prec_est()] and the post-processing step
#' [fact_prec_processing()] keep the interface of the original PrecFact code.
#'
#' @references
#' Chandra NK, Mueller P, Sarkar A (2026). Bayesian Scalable Precision Factor Analysis for
#' Gaussian Graphical Models. Bayesian Analysis, 21(1), 75-103.
#' \doi{10.1214/24-BA1461}
#'
#' @keywords internal
#' @useDynLib PrecFact, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

utils::globalVariables(c("name", "magnitude", "sign_label", "strength", "x", "y",
                         "Var1", "Var2", "value", "epsilon", "pFDR", "n_edges",
                         "iteration", "logdet", "chain"))
