#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# Usage:
# Rscript --vanilla plotTreeAndDomains3.r <treefile> <rps-blast.out> <homologs.fas> <output.pdf>

# --------------------
# Argument check
# --------------------
if (length(args) < 4) {
  stop("You must provide (1) a tree file (2) domain file (3) fasta sequences (4) output file name for pdf.", call. = FALSE)
}

# --------------------
# Libraries
# --------------------
suppressPackageStartupMessages({
  library(ggtree)
  library(data.table)
  library(drawProteins)
  library(ggplot2)
  library(seqinr)
  library(grid)  # for unit()
})

# Helper: define multiplot if missing (simple wrapper)
if (!exists("multiplot")) {
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, widths=NULL, heights=NULL, ncol=NULL, nrow=NULL) {
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
      if (!is.null(ncol)) cols <- ncol
      if (!is.null(nrow)) {
        layout <- matrix(seq(1, cols * nrow), nrow = nrow, ncol = cols, byrow = TRUE)
      } else {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, byrow = TRUE)
      }
    }
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),
                                                widths = if(is.null(widths)) unit(rep(1, ncol(layout)), "null") else widths,
                                                heights = if(is.null(heights)) unit(rep(1, nrow(layout)), "null") else heights)))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# --------------------
# Read and prep tree
# --------------------
hoxt <- read.tree(args[1])

# Tree plots (two variants, one with branch.length = 'none')
t  <- ggtree(hoxt, aes(x, y), branch.length = 'none') +
  geom_tree() + theme_tree() +
  geom_tiplab(cex = 2.5) +
  geom_label2(size = 2.5, aes(subset = !isTip, label = label),
              nudge_x = 0, nudge_y = .175, label.padding = unit(0.05, "lines")) +
  ggplot2::xlim(0, max(24, 0.389 * (sum(hoxt$edge.length)) - 26.389)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

t2 <- ggtree(hoxt, aes(x, y)) +
  geom_tree() + theme_tree() +
  geom_tiplab(cex = 2.5) +
  geom_label2(size = 2.5, aes(subset = !isTip, label = label),
              nudge_x = 0, nudge_y = .175, label.padding = unit(0.05, "lines")) +
  ggplot2::xlim(0, max(24, 0.389 * (sum(hoxt$edge.length)) - 26.389)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

torder <- data.frame("entryName" = get_taxa_name(t), "order" = length(get_taxa_name(t)):1)

# --------------------
# Read RPS-BLAST and FASTA, coerce numeric, build drawProteins table
# --------------------
# Expecting outfmt 6: qseqid qlen qstart qend evalue stitle (tab-separated)
hoxd <- fread(args[2], sep = "\t", header = FALSE, data.table = FALSE)

# Keep qseqid, qlen, qstart, qend, stitle
rel_data0 <- hoxd[, c(1:4, 6)]
colnames(rel_data0) <- c("V1","V2","V3","V4","V5")
rel_data0$type <- "DOMAIN"

# Coerce numeric where needed (qlen, qstart, qend)
suppressWarnings({
  rel_data0$V2 <- as.numeric(rel_data0$V2)  # qlen
  rel_data0$V3 <- as.numeric(rel_data0$V3)  # qstart
  rel_data0$V4 <- as.numeric(rel_data0$V4)  # qend
})

# Drop any malformed rows (non-numeric starts/ends)
rel_data0 <- rel_data0[!is.na(rel_data0$V3) & !is.na(rel_data0$V4), , drop = FALSE]

# CHAIN rows from fasta (full-length proteins)
reldatachain0 <- read.fasta(args[3])
reldatachain <- data.frame(
  V1 = names(reldatachain0),
  V2 = getLength(reldatachain0),
  V3 = 1,
  V4 = getLength(reldatachain0),
  V5 = names(reldatachain0),
  type = "CHAIN",
  stringsAsFactors = FALSE
)

# Combine DOMAIN + CHAIN
rel_data1 <- rbind(rel_data0, reldatachain)

# drawProteins expects: type, description, begin, end, length, accession, entryName, taxid
rel_data2 <- rel_data1[, c("type","V5","V3","V4","V2","V1","V1")]
colnames(rel_data2) <- c("type","description","begin","end","length","accession","entryName")
rel_data2$taxid <- 1L

# Ensure numeric after rbind
suppressWarnings({
  rel_data2$begin  <- as.numeric(rel_data2$begin)
  rel_data2$end    <- as.numeric(rel_data2$end)
  rel_data2$length <- as.numeric(rel_data2$length)
})

# Truncate long descriptions for clean legend
rel_data2$description <- substr(rel_data2$description, 1, 100)

# Keep only rows with valid coordinates (critical for draw_canvas)
rel_data2 <- rel_data2[!is.na(rel_data2$begin) & !is.na(rel_data2$end), , drop = FALSE]

# Merge plotting order
rel_data3 <- merge(rel_data2, torder, by = "entryName", all.x = TRUE)

# --------------------
# Guard against empty data (no domains/chain rows)
# --------------------
if (nrow(rel_data3) == 0) {
  warning("No valid domain/chain rows to draw; plotting tree only.")
  pdf(args[4], width = 10.5, height = 7)
  print(t)
  print(t2)
  dev.off()
  cat(paste0(args[4], ", a pdf file, has been outputted\n"))
  quit(save = "no")
}

# --------------------
# Draw proteins canvas + layers
# --------------------
p0 <- draw_canvas(rel_data3)
p1 <- draw_chains(p0, rel_data3, label_chains = FALSE, outline = "black")
p2 <- draw_domains(p1, rel_data3, label_domains = FALSE, show.legend = TRUE)

p <- draw_domains(p1, rel_data3, label_domains = FALSE, label_size = 0.7, show.legend = FALSE) +
  theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, 'cm'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(-8, 0, -10, -5), "mm"),
        legend.position = c(-0.25, 0.25),
        legend.background = element_rect(fill = "transparent"))

# --------------------
# Output PDF
# --------------------
pdf(args[4], paper = "USr", width = 10.5, height = 7)
multiplot(t,  p, ncol = 2, widths = c(8, 2))
multiplot(t2, p, ncol = 2, widths = c(8, 2))
p2
dev.off()

cat(paste0(args[4], ", a pdf file, has been outputted\n"))
