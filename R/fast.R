#' Fast inverse and log-determinant of a low-rank-plus-diagonal matrix
#'
#' For \eqn{\Omega = \Lambda\Lambda^T + \Delta} with \eqn{\Lambda} of size \eqn{d \times q},
#' \eqn{q \ll d}, and \eqn{\Delta = diag(ps)}, `fast_fact_inv()` returns \eqn{\Omega^{-1}}
#' by the Sherman-Morrison-Woodbury identity and `fast_fact_det()` returns
#' \eqn{\log\det\Omega}, both at a cost of order \eqn{d q^2} plus the \eqn{d^2} output.
#'
#' @usage
#' fast_fact_inv(Lambda, ps)
#' fast_fact_det(Lambda, ps)
#' @param Lambda A \eqn{d \times q} matrix.
#' @param ps A vector of \eqn{d} positive numbers, the diagonal of \eqn{\Delta}.
#'
#' @return `fast_fact_inv()`: a \eqn{d \times d} matrix. `fast_fact_det()`: a number.
#' @examples
#' set.seed(1)
#' L <- matrix(rnorm(12), 6, 2); ps <- runif(6, 1, 2)
#' Om <- tcrossprod(L) + diag(ps)
#' all.equal(fast_fact_inv(L, ps), solve(Om))
#' all.equal(fast_fact_det(L, ps), as.numeric(determinant(Om)$modulus))
#' @name fast_fact
#' @aliases fast_fact_inv fast_fact_det
#' @export fast_fact_inv
#' @export fast_fact_det
NULL
