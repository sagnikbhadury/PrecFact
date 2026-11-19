# README figures.
#   Rscript tools/make_readme_figures.R <dir with tnbc_fit.rds> <figure dir> [library]
args <- commandArgs(trailingOnly = TRUE)
res_dir <- args[1]; fig <- args[2]
lib <- if (length(args) > 2) args[3] else NULL
library(PrecFact, lib.loc = lib)
library(ggplot2)
library(patchwork)
save_png <- function(p, file, w, h) ggsave(file.path(fig, file), p, width = w, height = h, dpi = 150, bg = "white")
white <- theme(plot.background = element_rect(fill = "white", colour = NA))

# ---- Simulated example (design of the original vignette) ---------------------------
sim <- simulate_precfact(n = 200, d = 100, q = 8, pr = 0.1, val = 2.5, seed = 1)
fit_sim <- precfact(sim$y, chains = 2, seed = 1, verbose = FALSE)
print(fit_sim)
up <- upper.tri(sim$graph)
truth <- abs(sim$pcor[up]) > 0.1   # "true" edge if |partial correlation| > 0.1, as in the article
est <- fit_sim$graph[up]
cat(sprintf("Simulation: %d true edges (|pcor| > 0.1); %d selected; sensitivity %.3f, specificity %.4f, false discoveries %.3f; cor(est, true pcor) %.3f\n",
            sum(truth), sum(est), mean(est[truth]), mean(!est[!truth]),
            sum(est & !truth) / max(sum(est), 1), cor(fit_sim$pcor[up], sim$pcor[up])))
tm <- sim$pcor; diag(tm) <- NA; tm[abs(tm) < 1e-10] <- NA
em <- graph_pcor(fit_sim); diag(em) <- NA; em[em == 0] <- NA
lim <- max(abs(c(tm, em)), na.rm = TRUE)
sc <- scale_fill_gradient2(name = "Partial\ncorrelation", low = "#2166AC", mid = "white", high = "#C0392B",
                           midpoint = 0, limits = c(-lim, lim), na.value = "white")
p <- (heatplot(tm, show_labels = FALSE, diverging = TRUE) + sc + ggtitle("True")) +
  (heatplot(em, show_labels = FALSE, diverging = TRUE) + sc + ggtitle("PrecFact estimate (selected edges)")) +
  plot_layout(guides = "collect")
save_png(p, "sim_heatmaps.png", 9, 4.4)
tg <- sim$pcor; tg[abs(tg) < 0.1] <- 0
p <- (plot_network(tg, show_labels = FALSE, title = "True graph (|partial correlation| > 0.1)") |
        plot_network(fit_sim, show_labels = FALSE, title = "PrecFact graph")) &
  theme(legend.position = "none")
save_png(p, "sim_networks.png", 10, 5.2)

# ---- TNBC ---------------------------------------------------------------------------
f <- readRDS(file.path(res_dir, "tnbc_fit.rds"))
print(f)
p <- plot_network(f, layout = "stress", label_size = 2.1, node_size_range = c(0.5, 3.5),
                  edge_width_range = c(0.2, 1.6),
                  title = sprintf("METABRIC TNBC: PrecFact graph of %d genes", f$d),
                  subtitle = sprintf("%d edges among %d connected genes; 4 chains pooled; posterior FDR %.3f",
                                     sum(f$graph[upper.tri(f$graph)]), sum(rowSums(f$graph) > 0), f$pFDR))
save_png(p + theme(legend.box = "horizontal"), "tnbc_network_full.png", 9, 7.6)

hubs <- top_hubs(f)
min_deg <- 3
hub_genes <- hubs$node[hubs$degree >= min_deg]
cat("Hub genes (degree >= 3):", length(hub_genes), "\n")
p <- plot_network(f, nodes = hub_genes, label_size = 2.8,
                  title = "Hub genes of the TNBC graph",
                  subtitle = sprintf("Genes with at least %d edges, and the edges among them", min_deg))
save_png(p, "tnbc_network_hubs.png", 7.5, 8.2)

conn <- hubs$node
ord <- hclust(as.dist(1 - abs(graph_pcor(f, conn))), method = "average")$order
p <- heatplot(f, nodes = conn[ord], show_labels = TRUE) +
  ggtitle("Selected partial correlations among the connected genes") +
  theme(axis.text = element_text(size = 4.5))
save_png(p, "tnbc_heatmap.png", 9, 8.4)

top <- head(precfact_edges(f), 25)
top$pair <- factor(paste(top$from, "-", top$to), levels = rev(paste(top$from, "-", top$to)))
p <- ggplot(top, aes(x = partial_correlation, y = pair, colour = partial_correlation > 0)) +
  geom_errorbar(aes(xmin = partial_correlation - 2 * sd, xmax = partial_correlation + 2 * sd), width = 0,
                orientation = "y") +
  geom_point() + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2166AC"), guide = "none") +
  labs(x = "Partial correlation (posterior mean +/- 2 sd, pooled chains)", y = NULL,
       title = "The 25 strongest edges") + theme_minimal(base_size = 10) + white
save_png(p, "tnbc_top_edges.png", 7, 6)

ft <- f$fdr_table
p2 <- ggplot(ft, aes(epsilon, pFDR)) + geom_line() + geom_point(aes(size = n_edges)) +
  geom_hline(yintercept = f$settings$alpha, linetype = 2, colour = "#C0392B") +
  geom_vline(xintercept = f$epsilon, linetype = 3) +
  scale_size(name = "Edges", range = c(0.5, 3)) +
  labs(x = "epsilon (interval null |rho| <= epsilon)", y = "Posterior FDR",
       title = "Posterior FDR over the grid") + theme_minimal(base_size = 10)
p1 <- plot_trace(f) + ggtitle("Trace of log det(Omega), four chains")
save_png((p1 | p2) & white, "tnbc_diagnostics.png", 10, 3.8)
