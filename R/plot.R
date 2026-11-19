# Figures. These replace derive_graph() of the original code, which needed the
# 'colortools' package (no longer on CRAN), and keep heatplot().

.as_pcor <- function(x, type) {
  if (inherits(x, "precfact")) return(graph_pcor(x))
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x))
    stop("`x` must be a precfact object or a square numeric matrix.", call. = FALSE)
  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  rownames(x) <- colnames(x)
  if (type == "precision") {
    isd <- 1 / sqrt(diag(x))
    x <- -x * outer(isd, isd)
    diag(x) <- 1
  }
  x
}

#' Edge and node tables of a network
#'
#' @param x A `precfact` object, or a symmetric matrix of partial correlations or a
#'   precision matrix (see `type`).
#' @param type For a matrix `x`: `"partial_correlation"` or `"precision"`.
#' @param nodes Optional subset of variables (names or indices).
#' @param threshold Entries with smaller absolute value are not edges.
#' @param keep_isolated Keep nodes without edges?
#'
#' @return A list with `edges` (`from`, `to`, `partial_correlation`), `nodes` (`name`,
#'   `degree`, `strength`) and `graph`, an `igraph` object.
#' @examples
#' sim <- simulate_precfact(n = 50, d = 10, q = 2, seed = 1)
#' network_data(sim$Omega, type = "precision", threshold = 0.1)$edges
#' @export
network_data <- function(x, type = c("partial_correlation", "precision"), nodes = NULL,
                         threshold = 1e-8, keep_isolated = FALSE) {
  type <- match.arg(type)
  m <- .as_pcor(x, type)
  if (!is.null(nodes)) m <- m[nodes, nodes, drop = FALSE]
  diag(m) <- 0
  nm <- colnames(m)
  idx <- which(upper.tri(m) & abs(m) >= threshold, arr.ind = TRUE)
  edges <- data.frame(from = nm[idx[, 1]], to = nm[idx[, 2]], partial_correlation = m[idx],
                      stringsAsFactors = FALSE)
  act <- abs(m) >= threshold
  nd <- data.frame(name = nm, degree = rowSums(act), strength = rowSums(abs(m) * act),
                   stringsAsFactors = FALSE, row.names = NULL)
  if (!keep_isolated) nd <- nd[nd$degree > 0, , drop = FALSE]
  list(edges = edges, nodes = nd,
       graph = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nd))
}

#' Draw a network
#'
#' Draws the selected graph with `ggraph`. Edges are coloured by the sign of the posterior
#' mean partial correlation and their width follows its magnitude; nodes are sized by
#' strength, the sum of absolute partial correlations. Isolated variables are dropped
#' unless `keep_isolated = TRUE`.
#'
#' @inheritParams network_data
#' @param layout A `ggraph` layout: `"circle"` (default), `"stress"`, `"fr"`, `"kk"`, ...
#' @param title,subtitle Text of the figure.
#' @param show_labels Draw node labels? Use `FALSE` for large graphs.
#' @param label_size Label size.
#' @param node_colour Node colour.
#' @param edge_alpha Edge transparency.
#' @param edge_width_range,node_size_range Ranges of edge widths and node sizes.
#' @param curvature Edge curvature for the circle layout; 0 gives straight lines.
#' @param seed Seed for stochastic layouts.
#'
#' @return A `ggplot` object.
#' @examples
#' sim <- simulate_precfact(n = 50, d = 12, q = 2, seed = 1)
#' plot_network(sim$Omega, type = "precision", threshold = 0.05, title = "True graph")
#' @export
plot_network <- function(x, type = c("partial_correlation", "precision"), nodes = NULL,
                         threshold = 1e-8, keep_isolated = FALSE, layout = "circle",
                         title = NULL, subtitle = NULL, show_labels = TRUE, label_size = 2.6,
                         node_colour = "grey25", edge_alpha = 0.7,
                         edge_width_range = c(0.2, 1.5), node_size_range = c(0.6, 4.5),
                         curvature = 0.15, seed = 1) {
  type <- match.arg(type)
  net <- network_data(x, type, nodes, threshold, keep_isolated)
  base_theme <- ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom", legend.box = "vertical",
                   legend.title = ggplot2::element_text(size = 8, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  if (nrow(net$edges) == 0) {
    return(ggplot2::ggplot() + base_theme +
             ggplot2::labs(title = title, subtitle = paste(c(subtitle, "no edges"), collapse = ": ")))
  }
  g <- net$graph
  igraph::E(g)$magnitude <- abs(net$edges$partial_correlation)
  igraph::E(g)$sign_label <- factor(ifelse(net$edges$partial_correlation > 0, "Positive", "Negative"),
                                    levels = c("Positive", "Negative"))
  circle <- identical(layout, "circle")
  set.seed(seed)
  p <- ggraph::ggraph(g, layout = layout)
  aes_e <- ggplot2::aes(colour = sign_label, width = magnitude)
  p <- p + if (circle && curvature != 0) {
    ggraph::geom_edge_arc(aes_e, strength = curvature, alpha = edge_alpha)
  } else {
    ggraph::geom_edge_link(aes_e, alpha = edge_alpha)
  }
  p <- p +
    ggraph::scale_edge_width(name = "|Partial correlation|", range = edge_width_range) +
    ggraph::scale_edge_colour_manual(name = "Sign", values = c(Positive = "#C0392B", Negative = "#2166AC"),
                                     drop = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(size = strength), colour = node_colour) +
    ggplot2::scale_size(name = "Strength", range = node_size_range)
  if (show_labels) {
    if (circle) {
      p <- p + ggraph::geom_node_text(
        ggplot2::aes(label = name, x = x * 1.1, y = y * 1.1,
                     angle = ifelse(atan(-(x / y)) * (180 / pi) < 0,
                                    90 + atan(-(x / y)) * (180 / pi),
                                    270 + atan(-x / y) * (180 / pi)),
                     hjust = ifelse(x > 0, 0, 1)),
        size = label_size)
      p <- p + ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
    } else {
      p <- p + ggraph::geom_node_text(ggplot2::aes(label = name), size = label_size, repel = TRUE,
                                      segment.colour = "grey70", segment.size = 0.2)
    }
  }
  p + ggplot2::coord_fixed() + base_theme + ggplot2::labs(title = title, subtitle = subtitle)
}

#' Heat map of a matrix
#'
#' Draws a matrix as a heat map. For a `precfact` fit it shows the selected partial
#' correlation matrix ([graph_pcor()]), with zero entries left blank, on a diverging scale.
#' For a general matrix it keeps the behaviour of the original `heatplot()` (viridis scale,
#' `NA` entries blank).
#'
#' @param x A `precfact` object or a numeric matrix.
#' @param nodes Optional subset of variables for a fit.
#' @param show_labels Show axis labels?
#' @param diverging Use a blue-white-red scale centred at 0? Defaults to `TRUE` for a fit.
#'
#' @return A `ggplot` object.
#' @examples
#' sim <- simulate_precfact(n = 50, d = 20, q = 2, seed = 1)
#' C <- sim$Omega; C[C == 0] <- NA
#' heatplot(C)
#' @export
heatplot <- function(x, nodes = NULL, show_labels = NULL, diverging = inherits(x, "precfact")) {
  if (inherits(x, "precfact")) {
    m <- graph_pcor(x, nodes)
    diag(m) <- NA
    m[m == 0] <- NA
  } else {
    m <- as.matrix(x)
  }
  if (is.null(colnames(m))) colnames(m) <- seq_len(ncol(m))
  if (is.null(rownames(m))) rownames(m) <- seq_len(nrow(m))
  if (is.null(show_labels)) show_labels <- ncol(m) <= 60
  df <- data.frame(Var1 = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
                   Var2 = factor(rep(rownames(m), times = ncol(m)), levels = rev(rownames(m))),
                   value = as.vector(m))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var1, y = Var2, fill = value)) + ggplot2::geom_tile()
  p <- p + if (diverging) {
    lim <- max(abs(df$value), na.rm = TRUE)
    ggplot2::scale_fill_gradient2(name = "Partial\ncorrelation", low = "#2166AC", mid = "white",
                                  high = "#C0392B", midpoint = 0, limits = c(-lim, lim),
                                  na.value = "white")
  } else {
    ggplot2::scale_fill_viridis_c(option = "B", na.value = "white")
  }
  p <- p + ggplot2::coord_fixed() + ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!show_labels) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                           axis.text.y = ggplot2::element_blank())
  p
}

#' Trace plot of the log-determinant of the precision matrix
#'
#' A quick convergence check: the log-determinant of \eqn{\Omega} in each kept draw.
#'
#' @param fit A `precfact` object.
#' @return A `ggplot` object.
#' @examples
#' sim <- simulate_precfact(n = 80, d = 10, q = 2, seed = 1)
#' fit <- precfact(sim$y, nrun = 400, burn = 100, seed = 1, verbose = FALSE)
#' plot_trace(fit)
#' @export
plot_trace <- function(fit) {
  .check_fit(fit)
  ch <- if (is.null(fit$chain)) rep(1L, length(fit$logdet)) else fit$chain
  it <- stats::ave(seq_along(ch), ch, FUN = seq_along)
  df <- data.frame(iteration = fit$settings$burn + fit$settings$thin * it,
                   logdet = fit$logdet, chain = factor(ch))
  ggplot2::ggplot(df, ggplot2::aes(iteration, logdet, colour = chain)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = rep_len(c("grey25", "#E69F00", "#56B4E9", "#009E73",
                                                    "#CC79A7", "#D55E00"), nlevels(df$chain))) +
    ggplot2::labs(x = "MCMC iteration", y = expression(log ~ det ~ Omega), colour = "Chain") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}
