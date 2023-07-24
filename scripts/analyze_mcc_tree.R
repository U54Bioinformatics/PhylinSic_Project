library(argparse)

parser <- ArgumentParser()
parser$add_argument("filenames", nargs="+")
args <- parser$parse_args()

BEAST.FILE <- args$filenames[1]
METADATA.FILE <- args$filenames[2]
PHYLO.FILE <- args$filenames[3]
ROOTED.PHYLO.FILE <- args$filenames[4]
DIST.FILE <- args$filenames[5]
PHYLO.METADATA.FILE <- args$filenames[6]
ROOTED.PHYLO.METADATA.FILE <- args$filename[7]


library(phytools)
source("scripts/filelib.R")
source("scripts/beastlib.R")



# BEAST.FILE is combined MCC tree from Beast, NEXUS format.
# tidytree treedata S4 object.
mcc.tree <- treeio::read.beast(BEAST.FILE)

# Write out the cophenetic distances.
# cell x cell matrix of distances.  rownames and colnames are cell
# names.
cophenetic.dist <- cophenetic(mcc.tree@phylo)
x <- cophenetic.dist
data.out <- cbind(rownames(x), cophenetic.dist)
colnames(data.out) <- c("Cell", colnames(x))
my.write(data.out, DIST.FILE, col.names=TRUE)

# Read in the cell metadata and align to the cophenetic distances.
metadata <- my.read(METADATA.FILE)
if(!any(names(metadata) == "Cell"))
  stop('Metadata file does not have "Cell" column')
# May not have Outgroup.
#if(!any(names(metadata) == "Outgroup"))
#  stop('Metadata file does not have "Outgroup" column')
cells <- rownames(cophenetic.dist)
if(!all(cells == colnames(cophenetic.dist))) stop("unaligned")
I <- match(cells, metadata[["Cell"]])
missing <- cells[is.na(I)]
if(length(missing) > 0) {
  s <- ""
  if(length(missing) > 1) s <- "s"
  x <- sprintf("%d cell%s in tree not found in metadata.", length(missing), s)
  y <- missing
  if(length(y) <= 5) {
    y <- paste(y, collapse="\\n")
  } else {
    y <- paste(y[1:5], collapse="\\n")
    y <- sprintf("%s\\n[... +%d more]", y, length(missing)-5)
  }
  x <- sprintf("%s\\n%s", x, y)
  stop(x)
}
metadata <- metadata[I,]

# Re-root the tree, if there are any outgroup cells specified.

# Find the distances of the cells in the outgroup.
root.cell <- NULL
# If metadata doesn't have Outgroup, then "I" will be "integer(0)".
I <- which(metadata[["Outgroup"]] == "yes")
if(length(I) == 1){
  root.cell <- cells[I[1]]
} else if(length(I) >= 2) {
  # The root cell is the cell that has the minimum average distance to
  # the rest of the cells in the outgroup.
  outgroup.dist <- cophenetic.dist[I, I]
  m <- apply(outgroup.dist, 1, median)
  i <- which.min(m)
  root.cell <- cells[I[i]]
}

tree <- mcc.tree@phylo  # S3 object
rerooted.tree <- NULL
if(!is.null(root.cell)) {
  # resolve.root=TRUE will resolve as a bifurcating node.
  x <- ape::root(tree, outgroup=root.cell, resolve.root=TRUE)
  # No.  Preserve this.
  #MIN.DIST <- 1E-6
  #x$edge.length[x$edge.length < 0] <- MIN.DIST
  rerooted.tree <- x
}

# Write out the phylogeny in Newick format.
ape::write.tree(tree, PHYLO.FILE)
if(!is.null(rerooted.tree))
  ape::write.tree(rerooted.tree, ROOTED.PHYLO.FILE)

write.mcc.tree.metadata(mcc.tree, PHYLO.METADATA.FILE)
if(!is.null(rerooted.tree))
  write.tree.metadata(rerooted.tree, ROOTED.PHYLO.METADATA.FILE)
