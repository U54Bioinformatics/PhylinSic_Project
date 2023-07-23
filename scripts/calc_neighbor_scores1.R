REF.FILE <- snakemake@input[[1]]
ALT.FILE <- snakemake@input[[2]]
SCORE.FILE <- snakemake@output[[1]]
K <- snakemake@params$K
NUM.SAMPLES <- snakemake@params$num_samples
RNG.SEED <- snakemake@params$rng_seed
START <- snakemake@params$start
SKIP <- snakemake@params$skip
NUM.CORES <- snakemake@params$num_cores


library(parallel)
source("scripts/beastlib.R")


my.read <- function(filename, header=TRUE, skip=0, colClasses=NA, nrows=-1,
  clean.header=TRUE, blank.lines.skip=TRUE) {
  orig.filename <- filename
  if(length(grep("\\.gz$", filename, perl=TRUE, ignore.case=TRUE)))
    filename <- gzfile(orig.filename)
  D <- read.delim(filename, header=header, as.is=TRUE, comment.char="",
    quote="", skip=skip, colClasses=colClasses, nrows=nrows, 
    blank.lines.skip=blank.lines.skip)
  if(header & clean.header) {
    # Don't hash the header.
    filename <- orig.filename
    if(length(grep("\\.gz$", filename, perl=TRUE, ignore.case=TRUE)))
      filename <- gzfile(orig.filename)
    h <- read.delim(filename, header=FALSE, nrows=1, as.is=TRUE, skip=skip,
      comment.char="", quote="")
    names(D) <- h
  }
  D
}

my.write <- function(X, filename, row.names=FALSE, col.names=FALSE) {
  data.out <- as.matrix(X)
  if(is.logical(row.names)) {
    if(row.names)
      row.names <- rownames(X)
    else
      row.names <- c()
  }
  if(is.logical(col.names)) {
    if(col.names)
      col.names <- colnames(X)
    else
      col.names <- c()
  }
  if(length(col.names))
    data.out <- rbind(col.names, data.out)
  if(length(row.names)) {
    if(length(col.names))
      row.names <- c("", row.names)
    data.out <- cbind(row.names, data.out)
  }
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}




data.ref <- my.read(REF.FILE)
data.alt <- my.read(ALT.FILE)


if(nrow(data.ref) != nrow(data.alt)) stop("bad 1")
if(ncol(data.ref) != ncol(data.alt)) stop("bad 2")
if(any(data.ref[,1] != data.alt[,1])) stop("bad 3")
if(any(colnames(data.ref) != colnames(data.alt))) stop("bad 4")
  
M.ref <- as.matrix(data.ref[,2:ncol(data.ref)])
M.alt <- as.matrix(data.alt[,2:ncol(data.alt)])
M.ref[is.na(M.ref)] <- 0
M.alt[is.na(M.alt)] <- 0


# Calculate the probability distributions from cell j.
# m-length list of n x 3 matrices of probabilities.
p.j <- mclapply(1:ncol(M.ref), function(j) {
  calc.p.dist(M.ref[,j], M.alt[,j])
  }, mc.cores=NUM.CORES)

# This loop is really slow.
# Make a ncol x ncol matrix of cell similarities.  Optimization:
# calculate the upper diagonal only.

I <- seq(START, ncol(M.ref), SKIP)
x <- mclapply(I, function(j1) {
  set.seed(RNG.SEED)  # make sure this is reproducible.
  S.j <- rep(0, ncol(M.ref))
  # Will score lower diagonal of the matrix.
  for(j2 in j1:ncol(M.ref)) {
    S.j[j2] <- calc.cell.similarity(p.j[[j1]], p.j[[j2]], NUM.SAMPLES)
  }
  return(S.j)
}, mc.cores=NUM.CORES)
# Columns are the cells in this batch.
S <- matrix(unlist(x), nrow=ncol(M.ref), ncol=length(I))

# Write out the score file.
cell.names <- colnames(M.ref)
data.out <- cbind(cell.names, S)
colnames(data.out) <- c("Score", cell.names[I])
my.write(data.out, SCORE.FILE, col.names=TRUE)

