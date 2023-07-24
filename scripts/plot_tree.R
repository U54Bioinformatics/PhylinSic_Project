library(argparse)

parser <- ArgumentParser()
parser$add_argument("filenames", nargs="+")
args <- parser$parse_args()

TREE.FILE <- args$filenames[1]
METADATA.FILE <- args$filenames[2]
COLOR.BY.CATEGORY <- args$filenames[3]
HEIGHT <- as.numeric(args$filenames[4])
WIDTH <- as.numeric(args$filenames[5])
LAYOUT <- args$filenames[6]
OUTFILE <- args$filenames[7]

library(ape)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
source("scripts/filelib.R")
source("scripts/plotlib.R")

tree <- read.tree(TREE.FILE)
if(is.null(tree)) stop(sprintf("Problem reading tree file %s.", TREE.FILE))
metadata <- my.read(METADATA.FILE)

n <- max(tree$edge)
I <- match(tree$tip.label, metadata[["Cell"]])
if(any(is.na(I))) stop("tree metadata mismtch")
metadata <- metadata[I,]
category <- rep("", n)
category[1:length(I)] <- metadata[["Category"]]

# Really dumb and slow algorithm.
changed <- TRUE
while(changed) {
  changed <- FALSE
  for(i in 1:n) {
    I <- which(tree$edge[,1] == i)
    if(length(I) < 2)
      next
    j <- tree$edge[I,2]  # child nodes
    if(length(unique(category[j])) != 1)
      next
    if(category[i] == category[j][1])
      next
    category[i] <- category[j][1]
    changed <- TRUE
  }
}

p <- ggtree(tree, layout=LAYOUT)

if(COLOR.BY.CATEGORY) {
  all.categories <- sort(unique(category))
  palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(all.categories))
  x <- all.categories[2:length(all.categories)]
  p <- p + geom_tree(aes(color=category), layout=LAYOUT) +
    scale_color_manual(values=palette, limits=x)
  # Set the legend.
  p <- p + theme(legend.title=element_blank())  # no title
}

pdf(OUTFILE, height=HEIGHT, width=WIDTH)
print(p)
dev.off()
