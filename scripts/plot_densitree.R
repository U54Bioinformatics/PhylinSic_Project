library(argparse)

parser <- ArgumentParser()
parser$add_argument("filenames", nargs="+")
args <- parser$parse_args()

MODEL.FILE <- args$filenames[1]
I.START <- as.numeric(args$filenames[2])
I.END <- as.numeric(args$filenames[3])
HEIGHT <- as.numeric(args$filenames[4])
WIDTH <- as.numeric(args$filenames[5])
LAYOUT <- args$filenames[6]
OUTFILE <- args$filenames[7]

library(phangorn)
source("scripts/beastlib.R")

beast.model <- readRDS(MODEL.FILE)

if(I.START < 1) stop("bad 1")
if(I.START > length(beast.model$trees)) stop("bad 2")
if(I.END < 1) stop("bad 3")
if(I.END > length(beast.model$trees)) stop("bad 4")

# Make the consensus tree of these samples.
max.trees <- 100
trees <- beast.model$trees[I.START:I.END]
if(length(trees) > max.trees) {
  I <- seq(1, length(trees), length(trees)/max.trees)
  I <- round(I)
  I <- I[!duplicated(I)]
  trees <- trees[I]
}
if(0) {
  st <- superTree(trees, method="MRP", rooted=TRUE)
} else {
  # p is the proportion for a clade to be represented in the consensus
  # tree
  st <- consensus(trees, p=0.5)
}
x <- st$edge[,2]
x <- x[x <= length(st$tip.label)]
tips <- st$tip.label[x]

# Guess a good alpha.
# TREES  alpha
#  200    0.01
#   91    0.03    0.057 is too high.
#   40    0.06
#   20     0.1
#alpha <- 0.1-length(trees)/2300
alpha <- 1/length(trees)   # default
alpha <- min(max(alpha, 0.001), 1)
                                
#tree.type <- "cladogram"
tree.type <- "phylogram"

pdf(OUTFILE, height=HEIGHT, width=WIDTH, onefile=F)
par(mar=c(0, 0, 0, 0))
densiTree(trees, type=tree.type, width=1, scaleX=FALSE, cex=0.01,
  alpha=alpha, consensus=tips, scale.bar=FALSE)
dev.off()
