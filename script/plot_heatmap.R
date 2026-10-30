#!/usr/bin/env Rscript
# plot_heatmap.R
# Core SNP distance heatmap from snp-dists
# Author: Zhu Zhihao  Date: 2026.10.30
# Usage: Rscript plot_heatmap.R <snpdists.tsv> <outdir>
# Needs: ggplot2, scales

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript plot_heatmap.R <snpdists.tsv> <outdir>")
mat_path <- args[[1]]
outdir   <- args[[2]]
if (!file.exists(mat_path)) stop("Input not found: ", mat_path)

font_family <- if (.Platform$OS.type == "windows") "Arial" else "sans"
cluster_h <- 100
pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
         "#56B4E9", "#D55E00", "#F0E442", "#999999")

raw <- read.delim(mat_path, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
mat <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- as.character(raw[[1]])
colnames(mat) <- colnames(raw)[-1]

mat <- mat[rownames(mat) != "Reference", colnames(mat) != "Reference", drop = FALSE]
ids <- intersect(rownames(mat), colnames(mat))
if (length(ids) < 2L) stop("Need >= 2 samples after dropping Reference")
mat <- mat[ids, ids, drop = FALSE]

d <- (mat + t(mat)) / 2
diag(d) <- 0
hc <- hclust(as.dist(d), method = "average")
cl <- cutree(hc, h = cluster_h)
ord <- hc$labels[hc$order]
mat <- mat[ord, ord, drop = FALSE]
cl <- cl[ord]
n <- nrow(mat)

cid <- sort(unique(cl))
cl_lab <- setNames(LETTERS[seq_along(cid)], as.character(cid))
cl_col <- setNames(pal[((seq_along(cid) - 1L) %% length(pal)) + 1L], as.character(cid))

lab <- colnames(mat)
df <- expand.grid(
  SampleY = factor(lab, levels = rev(lab)),
  SampleX = factor(lab, levels = lab),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
df$SNPs <- as.vector(mat)
df_lab <- subset(df, SNPs <= cluster_h)

leg <- data.frame(
  SampleX = factor(lab[1], levels = lab),
  SampleY = factor(lab[1], levels = rev(lab)),
  Cluster = factor(unname(cl_lab[as.character(cid)]), levels = unname(cl_lab))
)

snp_max <- max(df$SNPs, na.rm = TRUE)
if (!is.finite(snp_max) || snp_max <= 0) snp_max <- 1
nice <- function(x) {
  if (x >= 1000) as.integer(round(x / 1000) * 1000) else as.integer(round(x))
}
br <- unique(c(0L, nice(snp_max / 2), nice(snp_max)))
fill_max <- max(snp_max, max(br))
axis_cex <- if (n <= 15) 9 else if (n <= 40) 7 else 5.5
cell_cex <- if (n <= 12) 3.2 else if (n <= 25) 2.4 else 1.8
w <- min(9, max(3.5, 0.32 * n + 2.0))
h <- min(8, max(3.2, 0.32 * n + 1.2))
x_cols <- unname(cl_col[as.character(cl[lab])])
y_cols <- unname(cl_col[as.character(cl[rev(lab)])])

suppressWarnings({
  p <- ggplot(df, aes(SampleX, SampleY)) +
    geom_tile(aes(fill = SNPs), color = "white", linewidth = 0.2) +
    geom_text(data = df_lab, aes(label = SNPs),
              colour = "black", size = cell_cex, family = font_family,
              show.legend = FALSE) +
    geom_point(data = leg, aes(colour = Cluster), size = 0, alpha = 0) +
    scale_fill_gradientn(
      name = "SNP",
      colours = c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5"),
      trans = "sqrt", breaks = br, labels = as.character(br),
      limits = c(0, fill_max),
      guide = guide_colorbar(
        barwidth = unit(0.35, "cm"), barheight = unit(3.2, "cm"),
        title.position = "top", title.hjust = 0.5, order = 1
      )
    ) +
    scale_colour_manual(
      name = "Cluster",
      values = setNames(unname(cl_col[as.character(cid)]), unname(cl_lab)),
      breaks = unname(cl_lab),
      guide = guide_legend(
        override.aes = list(size = 3.5, shape = 15, alpha = 1, stroke = 0),
        title.position = "top", title.hjust = 0, order = 2
      )
    ) +
    coord_fixed() +
    theme_bw(base_size = 8, base_family = font_family) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                size = axis_cex, colour = x_cols, face = "plain"),
      axis.text.y = element_text(size = axis_cex, colour = y_cols, face = "plain"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = axis_cex, face = "plain"),
      legend.text = element_text(size = axis_cex, face = "plain"),
      plot.margin = margin(4, 4, 4, 4)
    )
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
ggsave(file.path(outdir, "heatmap.pdf"), p, width = w, height = h, device = pdf_dev)
ggsave(file.path(outdir, "heatmap.png"), p, width = w, height = h, dpi = 300)
message("Wrote ", file.path(outdir, "heatmap.pdf"), " (", n, " samples, ", length(cid), " clusters)")
