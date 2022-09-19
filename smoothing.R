# FUNCTION:
# make.genotype.and.probability.matrices
#               Calls the genotypes across sites and cells.
#
# INPUTS:
# M.ref         A matrix of reference allele read counts.  Each row 
#               is a site, and each column is a cell.  Missing values 
#               should be filled with 0.
# M.alt         A matrix of alternate allele read counts.  Should be 
#               parallel to M.ref.
# K             The number of neighbors to use for smoothing.
#               Default is 10.
# DELTA         Amount of probability to assign to neighboring cells.
#               Should be between 0 (least smoothing) and 1 (most 
#               smoothing).  If 0, will do no smoothing.
# num.features  Number of sites to use when calculating neighbors.
#               By default, will use all sites.  If there are a lot 
#               of sites, will run a lot faster if you use fewer sites 
#               (e.g. 200).
# num.samples   Number of samples to use when comparing the genotype 
#               profiles of two cells.
#               Default is 100.
# sample        TRUE/FALSE.  If TRUE, will sample a genotype from its
#               probability distribution.  If FALSE (default), will select 
#               the highest probability genotype.
# rng.seed      Seed for random number generator.  Default is 0.
# num.cores     Number of cores to use.  Can specify more cores
#               to speed up computation.  Default is 1.
make.genotype.and.probability.matrices <- function(M.ref, M.alt,
  K=NULL, DELTA=NULL, num.features=NULL, num.samples=NULL, sample=FALSE,
  rng.seed=0, num.cores=1) {
  p.dist <- calc.p.dist.matrix(M.ref, M.alt, K=K, DELTA=DELTA,
    num.features=num.features, num.samples=num.samples, rng.seed=rng.seed,
    num.cores=num.cores)
  geno.mat <- matrix(NA, nrow(M.ref), ncol(M.ref))
  p.mat <- matrix(NA, nrow(M.ref), ncol(M.ref))
  for(i in 1:nrow(M.ref)) {
    p <- t(matrix(unlist(p.dist[[i]]), ncol=ncol(M.ref)))
    if(sample) {
      geno <- sample.genotype.from.p(p)
    } else {
      geno <- assign.genotype.from.p(p)
    }
    geno.mat[i,] <- geno
    p.mat[i,] <- apply(p, 1, max)
  }
  list(geno.mat, p.mat)
}

calc.p.dist.matrix <- function(M.ref, M.alt, K=NULL, DELTA=NULL,
  num.features=NULL, num.samples=NULL, rng.seed=0, num.cores=1) {
  # M.ref and M.alt are n x m matrices of the reference and alt
  # counts, where n is the number of sites, and m is the number of
  # cells.
  # Return an n x m x 3 list of lists of probability distributions.
  # K is the number of nearest neighbors to consider, and DELTA is the
  # amount of probability assigned to the neighboring cells.  0 does
  # no smoothing, and 1 does the most smoothing.
  # num.features is the number of features to use to calculate the
  # neighbors.
  # num.samples is the number of genotypes to sample when comparing
  # two cells.
  if(is.null(K))
    K <- 10
  if(is.null(DELTA))
    DELTA <- 0                      # by default, do no smoothing
  if(is.null(num.features))
    num.features <- nrow(M.ref)     # use all features by default
  if(is.null(num.samples))
    num.samples <- 100

  ## Choose the features with the highest VAF.
  #M.vaf <- M.alt / (M.ref + M.alt)
  #M.vaf[is.nan(M.vaf)] <- 0
  #x <- apply(M.vaf, 1, var)
  # Choose the features with the highest coverage.
  x <- apply(M.ref, 1, sum) + apply(M.alt, 1, sum)
  O <- order(x, decreasing=T)
  I <- O[1:num.features]

  # Calculate the probability distributions from cell j.
  # m-length list of n x 3 matrices of probabilities.
  p.j <- mclapply(1:ncol(M.ref), function(j) {
    calc.p.dist(M.ref[,j], M.alt[,j])
    }, mc.cores=num.cores)

  # If DELTA > 0, then smooth the probabilities using k nearest
  # neighbors.
  if(DELTA > 1E-6) {
    # Make a ncol x ncol matrix of cell similarities.  Optimization:
    # calculate the upper diagonal only.
    x <- mclapply(1:ncol(M.ref), function(j1) {
      set.seed(rng.seed)  # make sure this is reproducible.
      S.j <- rep(0, ncol(M.ref))
      for(j2 in j1:ncol(M.ref)) {
        S.j[j2] <- calc.cell.similarity(
          p.j[[j1]][I,], p.j[[j2]][I,], num.samples)
      }
      return(S.j)
    }, mc.cores=num.cores)
    S <- matrix(unlist(x), ncol(M.ref), ncol(M.ref))
    # Make symmetric.
    S <- pmax(S, t(S))

    # Find the nearest neighbors for each cell.
    NN <- matrix(NA, K, ncol(S))
    for(j in 1:ncol(S)) {
      O <- order(S[,j], decreasing=TRUE)
      NN[,j] <- sort(O[2:(2+K-1)])  # ignore highest score (matches itself)
    }

    # For each cell, average the probability distribution of the
    # neighbors.
    p.avgn <- list()   # m-length list of num.features x 3 matrices.
    for(j in 1:ncol(NN)) {
      p.avgn.j <- matrix(0, nrow(M.ref), 3)
      for(j2 in NN[,j])
        p.avgn.j <- p.avgn.j + p.j[[j2]]
      p.avgn[[j]] <- p.avgn.j / nrow(NN)
    }

    # Smooth the probability distribution.
    for(j in 1:length(p.j))
      p.j[[j]] <- (1-DELTA)*p.j[[j]] + DELTA*p.avgn[[j]]
  }

  dist.matrix <- list()
  for(i in 1:nrow(M.ref))
    dist.matrix[[i]] <- list()
  for(j in 1:length(p.j)) {
    p.jj <- p.j[[j]]
    colnames(p.jj) <- NULL
    for(i in 1:nrow(p.jj))
      dist.matrix[[i]][[j]] <- as.numeric(p.jj[i,])
  }
  dist.matrix
}

REF.CUT <- 0.30
ALT.CUT <- 0.70

calc.p.dist <- function(x.ref, x.alt) {
  # Calculate a probability distribution over HOMO-R, HET, HOMO-A.
  # for each cell.  x.ref and x.alt are n-length vectors of read
  # counts for each site.  Return an n x 3 matrix of probabilities.
  # Each row sums to 1.

  x1 <- pbeta(REF.CUT, x.alt+1, x.ref+1)
  x2 <- pbeta(ALT.CUT, x.alt+1, x.ref+1)
  p1 <- x1
  p2 <- x2 - x1
  p3 <- 1 - x2
  probs <- cbind(p1, p2, p3)

  probs
}

calc.cell.similarity <- function(p.geno.1, p.geno.2, num.samples) {
  # Return a score of the similarity between two cells.
  scores <- rep(NA, num.samples)
  for(i in 1:num.samples) {
    geno.1 <- sample.genotype.from.p(p.geno.1)
    geno.2 <- sample.genotype.from.p(p.geno.2)
    scores[i] <- sum(geno.1 == geno.2)
  }
  mean(scores) / nrow(p.geno.1) * 100
}

assign.genotype.from.p <- function(probs) {
  # Given a probability matrix from calc.p.dist, sample the genotypes.
  # Return a vector containing "R", "H", or "A" for homozygous
  # reference, heterozygous, or homozygous alt.
  x <- apply(probs, 1, which.max)
  geno <- rep(NA, nrow(probs))
  geno[x == 1] <- "R"
  geno[x == 2] <- "H"
  geno[x == 3] <- "A"
  geno
}

sample.genotype.from.p <- function(probs) {
  # Given a probability matrix from calc.p.dist, sample the genotypes.
  # Return a vector containing "R", "H", or "A" for homozygous
  # reference, heterozygous, or homozygous alt.
  p <- runif(nrow(probs))
  geno <- rep("A", nrow(probs))
  geno[p<(probs[,1]+probs[,2])] <- "H"
  geno[p<probs[,1]] <- "R"
  geno
}

make.dropout.matrix <- function(M.ref, M.alt) {
  # Return a matrix of 1/0, indicating whether that position is
  # dropped out.
  m <- matrix(0, nrow(M.ref), ncol(M.alt))
  m[(M.ref == 0) & (M.alt == 0)] <- 1
  m
}

