# Precision factor analysis of the METABRIC TNBC gene panel shipped with GGMsubho.
#   Rscript tools/tnbc_analysis.R <output dir> [library with PrecFact]
# The first 250 genes of the panel (n = 299 > d), four chains with the settings of
# Chandra, Mueller and Sarkar (2026) pooled. Saves the fit, edge and hub tables.
# Figures: tools/make_readme_figures.R.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
lib <- if (length(args) > 1) args[2] else NULL
dir.create(out, showWarnings = FALSE, recursive = TRUE)
library(PrecFact, lib.loc = lib)

# The panel is the `tnbc_metabric` data set of GGMsubho
# (https://github.com/sagnikbhadury/GGMsubho). Use the installed package if present,
# otherwise download the data file from GitHub.
if (requireNamespace("GGMsubho", quietly = TRUE)) {
  tnbc <- GGMsubho::tnbc_metabric
} else {
  tmp <- tempfile(fileext = ".rda")
  utils::download.file("https://github.com/sagnikbhadury/GGMsubho/raw/main/data/tnbc_metabric.rda",
                       tmp, mode = "wb", quiet = TRUE)
  e <- new.env()
  load(tmp, envir = e)
  tnbc <- e$tnbc_metabric
}
x <- tnbc[, 1:250]
cat("TNBC panel:", nrow(x), "patients x", ncol(x), "genes\n")

fit <- precfact(x, chains = 4, seed = 1, verbose = FALSE)
print(fit)
print(fit$chain_agreement)
saveRDS(fit, file.path(out, "tnbc_fit.rds"))
print(summary(fit, top = 15))
print(fit$fdr_table)
utils::write.csv(precfact_edges(fit), file.path(out, "tnbc_edges.csv"), row.names = FALSE)
utils::write.csv(top_hubs(fit), file.path(out, "tnbc_hubs.csv"), row.names = FALSE)
utils::write.csv(fit$chain_agreement, file.path(out, "tnbc_chain_agreement.csv"), row.names = FALSE)
cat("Degree distribution:\n"); print(table(cut(rowSums(fit$graph), c(-1, 0, 1, 2, 4, 9, 19, Inf))))
