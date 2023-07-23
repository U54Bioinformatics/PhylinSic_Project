REF.FILE <- snakemake@input[[1]]
ALT.FILE <- snakemake@input[[2]]
NEIGHBOR.FILE <- snakemake@input[[3]]

GENOTYPE.FILE <- snakemake@output[[1]]
PROBABILITY.FILE <- snakemake@output[[2]]

DELTA <- snakemake@params$delta
SAMPLE <- FALSE
IMPUTE <- snakemake@params$impute
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






if(is.null(DELTA))
  DELTA <- 0                      # by default, do no smoothing

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

# n           num mutations
# m           num cells
# K           number of neighbors
# M.ref       nxm matrix of counts
# M.alt       nxm matrix of counts
# M.neighbor  mxk matrix of indexes of neighbors
# p.j         m-length list of nx3 matrix
# p.avgn      m-length list of nx3 matrix

# Read the neighbors and convert into indexes.
data.neighbor <- my.read(NEIGHBOR.FILE)
cell.names <- colnames(M.ref)
M.neighbor <- matrix(NA, nrow=nrow(data.neighbor), ncol=ncol(data.neighbor))
for(i in 1:ncol(data.neighbor)) {
  I <- match(data.neighbor[,i], cell.names)
  if(any(is.na(I))) stop("bad")
  M.neighbor[,i] <- I
}



# Calculate the probability distributions from cell j.
p.j <- mclapply(1:ncol(M.ref), function(j) {
  calc.p.dist(M.ref[,j], M.alt[,j])
  }, mc.cores=NUM.CORES)

if(DELTA > 0) {
  # For each cell, average the probability distribution of the
  # neighbors.
  p.avgn <- list()
  for(j in 1:nrow(M.neighbor)) {
    p.avgn.j <- matrix(0, nrow(M.ref), 3)
    for(j2 in M.neighbor[j, 2:ncol(M.neighbor)])
      p.avgn.j <- p.avgn.j + p.j[[j2]]
    p.avgn[[j]] <- p.avgn.j / (ncol(M.neighbor)-1)
  }

  # Smooth the probability distribution.
  for(j in 1:length(p.j))
    p.j[[j]] <- (1-DELTA)*p.j[[j]] + DELTA*p.avgn[[j]]
}


M.geno <- matrix(NA, nrow(M.ref), ncol(M.ref))
M.p <- matrix(NA, nrow(M.ref), ncol(M.ref))
for(j in 1:length(p.j)) {
  p <- p.j[[j]]
  if(SAMPLE) {
    geno <- sample.genotype.from.p(p)
  } else {
    geno <- assign.genotype.from.p(p)
  }
  M.geno[,j] <- geno
  M.p[,j] <- apply(p, 1, max)
}

# Convert from fractions to probabilities.
M.p <- M.p * 100

if(!IMPUTE) {
  M.dropout <- make.dropout.matrix(M.ref, M.alt)
  M.geno[M.dropout==1] <- NA
}

# Write out the matrices.
x1 <- cbind(data.ref[,1], M.geno)
x2 <- cbind(data.ref[,1], M.p)
x1 <- rbind(colnames(data.ref), x1)
x2 <- rbind(colnames(data.ref), x2)
x1[is.na(x1)] <- ""
if(any(is.na(x2))) stop("bad")
write.table(x1, GENOTYPE.FILE, quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)
write.table(x2, PROBABILITY.FILE, quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)
