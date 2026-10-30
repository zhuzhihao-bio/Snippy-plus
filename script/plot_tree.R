#!/usr/bin/env Rscript
# plot_tree.R — core SNP tree (cluster color + support stars)
# Author: Zhu Zhihao  Date: 2026.10.30
# Usage: Rscript plot_tree.R <tree> <snpdists.tsv> <outdir>
# Cluster cut = heatmap (SNP <= 100); star if support >= 0.7 (or 70)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(ape)
})
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L)
  stop("Usage: Rscript plot_tree.R <tree> <snpdists.tsv> <outdir>")
tree_path <- args[[1]]
mat_path  <- args[[2]]
outdir    <- args[[3]]
if (!file.exists(tree_path)) stop("Tree not found: ", tree_path)
if (!file.exists(mat_path))  stop("Matrix not found: ", mat_path)

font_family <- if (.Platform$OS.type == "windows") "Arial" else "sans"
cluster_h <- 100
boot_cut  <- 0.7
pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
         "#56B4E9", "#D55E00", "#F0E442", "#999999")

# --- cluster from snp-dists ---
raw <- read.delim(mat_path, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
mat <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- as.character(raw[[1]])
colnames(mat) <- colnames(raw)[-1]
mat <- mat[rownames(mat) != "Reference", colnames(mat) != "Reference", drop = FALSE]
ids <- intersect(rownames(mat), colnames(mat))
mat <- mat[ids, ids, drop = FALSE]

d <- (mat + t(mat)) / 2
diag(d) <- 0
hc <- hclust(as.dist(d), method = "average")
cl <- cutree(hc, h = cluster_h)
cid <- sort(unique(cl))
cl_lab <- setNames(LETTERS[seq_along(cid)], as.character(cid))
cl_col <- setNames(pal[((seq_along(cid) - 1L) %% length(pal)) + 1L], as.character(cid))

# --- tree ---
tr <- read.tree(tree_path)
drop <- intersect(tr$tip.label, c("Reference", "reference"))
if (length(drop)) tr <- drop.tip(tr, drop)

tip_df <- data.frame(
  label = names(cl),
  Cluster = unname(cl_lab[as.character(cl)]),
  stringsAsFactors = FALSE
)
tip_df <- tip_df[tip_df$label %in% tr$tip.label, , drop = FALSE]
tip_df$Cluster <- factor(tip_df$Cluster, levels = unname(cl_lab))

boot_min <- boot_cut
if (!is.null(tr$node.label)) {
  supp <- suppressWarnings(as.numeric(tr$node.label))
  if (any(is.finite(supp) & supp > 1.5)) boot_min <- boot_cut * 100
}

n <- length(tr$tip.label)
cex <- if (n <= 15) 5 else if (n <= 40) 4 else 3.2
w <- min(9, max(5.5, 0.22 * n + 5))
h <- min(9, max(4, 0.32 * n + 2))

p <- ggtree(tr, size = 0.5) %<+% tip_df
xmax <- max(p$data$x[p$data$isTip], na.rm = TRUE)
tips <- p$data[p$data$isTip, , drop = FALSE]

p <- p +
  geom_segment(data = tips, aes(x = x, xend = xmax, y = y, yend = y),
               linetype = "dotted", linewidth = 0.3, colour = "grey60",
               inherit.aes = FALSE) +
  geom_text(data = tips, aes(x = xmax, y = y, label = label, colour = Cluster),
            hjust = 0, nudge_x = xmax * 0.02, size = cex,
            family = font_family, inherit.aes = FALSE, show.legend = FALSE) +
  geom_text2(
    aes(subset = !isTip & !is.na(as.numeric(label)) &
                  as.numeric(label) >= boot_min,
        label = "★"),
    colour = "red", size = 5, vjust = 0.35, family = font_family
  ) +
  scale_colour_manual(
    values = setNames(unname(cl_col[as.character(cid)]), unname(cl_lab)),
    guide = "none"
  ) +
  geom_treescale(x = 0, y = n + 0.6, width = signif(xmax * 0.15, 2),
                 fontsize = 4, family = font_family, offset = -0.4) +
  theme(plot.margin = margin(10, 20, 10, 10)) +
  xlim(0, xmax * 1.45)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
ggsave(file.path(outdir, "tree.pdf"), p, width = w, height = h, device = pdf_dev)
ggsave(file.path(outdir, "tree.png"), p, width = w, height = h, dpi = 300)
