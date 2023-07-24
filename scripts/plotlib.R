# Functions:
# green.shade
# red.shade
# white.shade
# rg.array.colors       Red/green array.
# by.array.colors       Blue/yellow array.
# red.green.soft
# red.blue.soft
# rgb.colors            From red, green, blue colorwheel.
# ryb.colors            From red, yellow, blue colorwheel.
# matlab.colors         Default Matlab "jet" colors.
# matlab.hot.colors
# broad.colors
# bild.colors           Andrea's version of the "jet" colors.
# yahoo.weather.colors  Like Yahoo Weather Map.
# genespring.colors
# 
# brewer.prgn.div
# brewer.rdbu.div
# brewer.rdylbu.div
# brewer.rdylgn.div
# brewer.spectral.div
# brewer.greens.seq
# brewer.blues.seq
# brewer.qual.set1
# brewer.qual.set2
# brewer.qual.set3
#
# get.palette.fn
# 
# my.lineplot
# my.lineplot2
# my.heatmap
# my.colorbar
# sample.colorbars      Print out all different color schemes.
# my.groupplot
# my.pcaplot
#
# my.violin
# my.beeswarm
# my.boxplot
# my.barplot
#
# get.bitmap.type
# sciNotation
# calc.plot.margins
#
# INTERNAL FUNCTIONS
# normalize.one.mv
# normalize.mv
# .groupplot.calc.X
# matrix2color
# matrix2rgb



BREWER_S1_RED <- "#D23329"
BREWER_S1_BLUE <- "#467CB4"
BREWER_S1_GREEN <- "#69AD59"
BREWER_S1_PURPLE <- "#8F559D"
BREWER_S1_ORANGE <- "#E97F32"
BREWER_S1_YELLOW <- "#FBFD60"
BREWER_S1_BROWN <- "#975C2C"
BREWER_S1_PINK <- "#E383BA"


green.shade <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 0),
    c(1.0,   0, 255, 0)),
    4, 2))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

red.shade <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0, 0, 0),
    c(1.0, 255, 0, 0)),
    4, 2))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

white.shade <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0, 0, 0),
    c(1.0, 255, 255, 255)),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

rg.array.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0, 255, 0),
    c(0.5,   0,   0, 0),
    c(1.0, 255,   0, 0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

by.array.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 255),
    c(0.5,   0,   0,   0),
    c(1.0, 255, 255,   0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

red.green.soft <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,   52, 112,  11),
    c(0.25,   95, 173,  93),
    c(0.50,  245, 245, 245),
    #c(0.75,  229,  89,  95),
    #c(1.00,  249,  39,  39),
    c(0.75,  229,  115,  105),
    c(1.00,  128,  52,  29)),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

red.blue.soft <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,    8,  34,  79),
    c(0.25,   79, 148, 194),
    c(0.50,  245, 245, 245),
    c(0.75,  213,  96,  76),
    c(1.00,   86,  15,  24)),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

rgb.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,   0,   0, 255),
    c(0.25,   0, 255, 255),
    c(0.50,   0, 255,   0),
    c(0.75, 255, 255,   0),
    c(1.00, 255,   0,   0)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

ryb.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,   0,   0, 255),
    c(0.25,   0, 255,   0),
    c(0.50, 255, 255,   0),
    c(0.75, 255, 128,   0),
    c(1.00, 255,   0,   0)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

matlab.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,   0,   0, 143),
    c(0.125,   0,   0, 255),
    c(0.250,   0, 127, 255),
    c(0.375,   0, 255, 255),
    c(0.500, 127, 255, 127),
    c(0.625, 255, 255,   0),
    c(0.750, 255, 127,   0),
    c(0.875, 255,   0,   0),
    c(1.000, 127,   0,   0)),
    4, 9))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

matlab.hot.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  10,   0,   0),
    c(0.361, 255,   0,   0),
    c(0.377, 255,  10,   0),
    c(0.738, 255, 255,   0),
    c(1.000, 255, 255, 255)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

broad.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  69,   0, 173),
    c(0.091,  39,   0, 209),
    c(0.182, 107,  88, 239),
    c(0.273, 136, 136, 255),
    c(0.364, 199, 193, 255),
    c(0.455, 213, 213, 255),
    c(0.545, 255, 192, 229),
    c(0.636, 255, 137, 137),
    c(0.727, 255, 112, 128),
    c(0.818, 255,  90,  90),
    c(0.909, 239,  64,  64),
    c(1.000, 214,  12,   0)),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

bild.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  49,  50, 114),
    c(0.050,  61,  69, 137),
    c(0.100,  62,  84, 154),
    c(0.150,  67,  89, 160),
    c(0.200,  85, 108, 176),
    c(0.250, 115, 145, 201),
    c(0.300, 160, 205, 240),
    c(0.350, 180, 220, 243),
    c(0.400, 169, 216, 211),
    c(0.450, 160, 208, 164),
    c(0.500, 179, 213, 112),
    c(0.550, 203, 220,  61),
    c(0.600, 232, 231,  61),
    c(0.650, 255, 234,  47),
    c(0.700, 250, 180,  50),
    c(0.750, 243, 136,  54),
    c(0.800, 231,  80,  61),
    c(0.850, 218,  54,  55),
    c(0.900, 204,  55,  59),
    c(0.950, 160,  52,  52),
    c(1.000, 114,  39,  44)),
    4, 21))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

yahoo.weather.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0, 255, 255, 255),
    c(0.1, 204, 255, 255),
    c(0.2, 153, 255, 255),
    c(0.3, 102, 204, 255),
    c(0.4,  84, 169, 255),
    c(0.5, 204, 255, 103),
    c(0.6, 255, 255, 103),
    c(0.7, 255, 204, 102),
    c(0.8, 255, 153, 102),
    c(0.9, 204, 102, 102),
    c(1.0, 209,  73,  73)),
    4, 11))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

genespring.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 255),
    c(0.5, 255, 255,   0),
    c(1.0, 255,   0,   0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.prgn.div <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 0, 69, 30),
    c(0.10, 26, 128, 63),
    c(0.20, 79, 174, 106),
    c(0.30, 169, 221, 162),
    c(0.40, 217, 241, 211),
    c(0.50, 248, 248, 248),
    c(0.60, 232, 212, 233),
    c(0.70, 195, 166, 208),
    c(0.80, 154, 112, 171),
    c(0.90, 119, 44, 135),
    c(1.00, 61, 1, 80)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.rdbu.div <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 6, 46, 95),
    c(0.10, 33, 102, 173),
    c(0.20, 67, 147, 197),
    c(0.30, 146, 198, 223),
    c(0.40, 210, 230, 240),
    c(0.50, 248, 248, 248),
    c(0.60, 253, 220, 200),
    c(0.70, 245, 165, 130),
    c(0.80, 215, 96, 78),
    c(0.90, 177, 20, 42),
    c(1.00, 103, 0, 28)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.rdylbu.div <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 44, 52, 147),
    c(0.10, 69, 118, 181),
    c(0.20, 116, 174, 209),
    c(0.30, 171, 217, 235),
    c(0.40, 227, 245, 249),
    c(0.50, 254, 255, 182),
    c(0.60, 254, 226, 144),
    c(0.70, 254, 175, 100),
    c(0.80, 245, 107, 70),
    c(0.90, 216, 47, 39),
    c(1.00, 154, 1, 64)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.rdylgn.div <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 10, 111, 54),
    c(0.10, 26, 152, 76),
    c(0.20, 103, 190, 99),
    c(0.30, 164, 217, 106),
    c(0.40, 218, 240, 140),
    c(0.50, 254, 255, 192),
    c(0.60, 254, 225, 139),
    c(0.70, 254, 175, 100),
    c(0.80, 245, 107, 70),
    c(0.90, 216, 47, 39),
    c(1.00, 157, 1, 55)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.spectral.div <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 102, 81, 164),
    c(0.10, 61, 143, 194),
    c(0.20, 102, 195, 165),
    c(0.30, 169, 221, 163),
    c(0.40, 231, 246, 152),
    c(0.50, 254, 255, 192),
    c(0.60, 254, 225, 139),
    c(0.70, 254, 175, 100),
    c(0.80, 245, 107, 70),
    c(0.90, 213, 62, 80),
    c(1.00, 150, 0, 69)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.reds.seq <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 255, 246, 242),
    c(0.14, 244, 223, 212),
    c(0.29, 240, 189, 166),
    c(0.43, 236, 145, 109),
    c(0.57, 229, 117, 80),
    c(0.71, 216, 79, 56),
    c(0.86, 187, 43, 34),
    c(1.00, 137, 27, 21)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.greens.seq <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 250, 253, 254),
    c(0.14, 232, 244, 249),
    c(0.29, 210, 235, 230),
    c(0.43, 166, 214, 200),
    c(0.57, 123, 191, 166),
    c(0.71, 98, 173, 126),
    c(0.86, 70, 138, 78),
    c(1.00, 35, 84, 39)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

brewer.blues.seq <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00, 248, 251, 252),
    c(0.14, 225, 234, 244),
    c(0.29, 203, 216, 234),
    c(0.43, 165, 199, 220),
    c(0.57, 125, 172, 208),
    c(0.71, 90, 144, 193),
    c(0.86, 71, 114, 176),
    c(1.00, 30, 67, 144)
    ),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

BREWER.S1.RED <- "#D23329"
BREWER.S1.BLUE <- "#467CB4"
BREWER.S1.GREEN <- "#69AD59"
BREWER.S1.PURPLE <- "#8F559D"
BREWER.S1.ORANGE <- "#E97F32"
BREWER.S1.YELLOW <- "#FBFD60"
BREWER.S1.BROWN <- "#975C2C"
BREWER.S1.PINK <- "#E383BA"

brewer.qual.set1 <- function(n) {
  COLORS <- c(
    BREWER.S1.RED,
    BREWER.S1.BLUE,
    BREWER.S1.GREEN,
    BREWER.S1.PURPLE,
    BREWER.S1.ORANGE,
    BREWER.S1.YELLOW,
    BREWER.S1.BROWN,
    BREWER.S1.PINK
  )
  if(n <= 0) stop("n must be greater than 0")
  if(n > length(COLORS)) stop(sprintf("n must be <= %d", length(COLORS)))
  return(COLORS[1:n])
}

brewer.qual.set2 <- function(n) {
  COLORS <- c(
    "#69C1A5",
    "#FB8C65",
    "#8DA1CA",
    "#EA89C1",
    "#A8D65B",
    "#FFD63F",
    "#E5C396",
    "#B0B0B0"
  )
  if(n <= 0) stop("n must be greater than 0")
  if(n > length(COLORS)) stop(sprintf("n must be <= %d", length(COLORS)))
  return(COLORS[1:n])
}

brewer.qual.set3 <- function(n) {
  COLORS <- c(
    "#8FD3C7",
    "#FFFDB6",
    "#BDBCDA",
    "#FA7F74",
    "#80B2D2",
    "#FCB267",
    "#B5DC6E",
    "#F7CCE4"
  )
  if(n <= 0) stop("n must be greater than 0")
  if(n > length(COLORS)) stop(sprintf("n must be <= %d", length(COLORS)))
  return(COLORS[1:n])
}


get.palette.fn <- function(name) {
  if(name == "red")
    return(red.shade)
  else if(name == "white")
    return(white.shade)
  else if(name == "red-green")
    return(rg.array.colors)
  else if(name == "blue-yellow")
    return(by.array.colors)
  else if(name == "matlab")
    return(matlab.colors)
  else if(name == "bild")
    return(bild.colors)
  else if(name == "genepattern")
    return(broad.colors)
  else if(name == "genespring")
    return(genespring.colors)
  else if(name == "yahoo")
    return(yahoo.weather.colors)
  else if(name == "brewer-prgn-div")
    return(brewer.prgn.div)
  else if(name == "brewer-rdbu-div")
    return(brewer.rdbu.div)
  else if(name == "brewer-rdylbu-div")
    return(brewer.rdylbu.div)
  else if(name == "brewer-rdylgn-div")
    return(brewer.rdylgn.div)
  else if(name == "brewer-spectral-div")
    return(brewer.spectral.div)
  stop(sprintf("Unknown color scheme: %s", name))
}


my.lineplot <- function(coords, xlim=NA, ylim=NA, col=NA, lwd=2) {
  # coords is a matrix where the columns are X_1, Y_1, ... X_N, Y_N.
  X <- coords[,seq(1, ncol(coords), 2)]
  Y <- coords[,seq(2, ncol(coords), 2)]
  if(ncol(X) != ncol(Y)) stop("X and Y unaligned")
  if(is.na(xlim))
    xlim <- c(min(X), max(X))
  if(is.na(ylim))
    ylim <- c(min(Y), max(Y))
  if(is.na(col))
    col <- matlab.colors(ncol(X))
  if(length(col) != ncol(X)) stop("colors unaligned")

  plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  for(i in 1:ncol(X)) {
    x <- X[,i]; y <- Y[,i]
    lines(x, y, lwd=lwd, col=col[i])
    points(x, y, col=col[i], pch=".", cex=8)
  }
}

my.lineplot2 <- function(X, Y=NA, xlab=NA, ylab=NA, 
  col=NA, SPACING=NA, lwd=1, main=NA) {
  # X should be gene x samples matrix.
  # Y should be a vector of the line number for each row of X (1-based).
  # col is a vector of the colors for each gene.
  # By default, prints from the bottom (row 1) up.
  if((length(Y) == 1) && is.na(Y))
    Y <- 1:nrow(X)
  num.lines <- max(Y)
  if((length(col)==1) && is.na(col))
    col <- rep("#000000", nrow(X))
  if((length(ylab)==1) && is.na(ylab))
    ylab <- NA
  if((length(xlab)==1) && is.na(xlab))
    xlab <- NA
  if((length(SPACING) == 1) && is.na(SPACING))
    SPACING <- max(abs(X))

  xlim <- c(1, ncol(X))  
  ylim <- c(-SPACING, num.lines*SPACING)
  plot(NA, type="n", axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="", 
    main=main)
  for(i in 1:nrow(X)) {
    offset <- (Y[i]-1) * SPACING
    x <- as.numeric(X[i,]) + offset
    lines(x, lwd=lwd, col=col[i])
  }
  axis(1, at=1:ncol(X), labels=xlab)
  axis(2, at=seq(0, (num.lines-1)*SPACING, SPACING), labels=ylab)
  box()
}

# matrix should contain values from 0 to 1.
my.heatmap <- function(matrix, col=rg.array.colors, xlab="", ylab="", 
  normalize=FALSE, scale=FALSE, cluster=FALSE) {
  # If normalize is TRUE, then will normalize to N(0, 1).  It can also
  # be a vector of (mean, variance) to specify the normalization
  # parameters.
  # If scale is TRUE, then will everything onto a 0 to 1 scale.  It
  # can also be a vector of (min, max, median) that indicates the
  # minimum, maximum, and median values that should correspond to 0,
  # 1, and 0.5.  If median is NA, then will not attempt to center the
  # median.
  # If cluster is TRUE, will cluster the rows and columns.
  if((length(normalize) == 1 && normalize == TRUE) || 
    (length(normalize) == 2)) {
    M <- 0; V <- 1
    if(length(normalize) == 2) {
      M <- normalize[1]; V <- normalize[2] 
    }
    m <- apply(matrix, 1, mean) - M
    matrix <- sweep(matrix, 1, m)
    matrix <- t(apply(matrix, 1, function(x) {
      V.0 <- var(x); M.0 <- mean(x); (x-M.0)*sqrt(V/V.0) + M.0 }))
  }

  if((length(scale) == 1 && scale == TRUE) || (length(scale) == 3)) {
    scale.min <- NA; scale.max <- NA; scale.med <- NA
    if(length(scale) != 1) {
      scale.min <- scale[1]; scale.max <- scale[2]; scale.med <- scale[3]
    }
    if(is.na(scale.min)) scale.min <- min(matrix)
    if(is.na(scale.max)) scale.max <- max(matrix)
    if(scale.max <= scale.min) stop("invalid scale parameters")
    matrix[matrix < scale.min] <- scale.min
    matrix[matrix > scale.max] <- scale.max
    if(!is.na(scale.med)) {
      # Center around 0, then scale so it's from -0.5 to +0.5.
      matrix <- matrix - scale.med
      x <- max(abs(c(min(matrix), max(matrix))))
      matrix <- matrix / x / 2 + 0.5
    } else {
      matrix <- matrix - min(matrix)
      matrix <- matrix / max(matrix)
    }
  }

  # col should be dark to bright.
  if(mode(col) == "character") {
    num.colors <- length(col)
  } else if(mode(col) == "list") {
    num.colors <- 256
    col <- rg.array.colors(num.colors)
  } else {
    num.colors <- 256
    col <- col(num.colors)
  }
  x.size <- 1
  y.size <- 1
  matrix <- as.matrix(matrix)

  # image treats the rows as x and columns as y.  So I need to
  # "rotate" the matrix 90 degrees clockwise.
  matrix <- matrix(t(matrix)[,nrow(matrix):1], ncol(matrix), nrow(matrix))

  # Weird.  If X11 is not running, then have to do it this way.
  # image puts row 1 on the bottom and column 1 on the right.  Flip
  # both the rows and columns.
  #matrix <- matrix[nrow(matrix):1,ncol(matrix):1]

  if(cluster) {
    h.r <- hclust(dist(matrix))
    h.c <- hclust(dist(t(matrix)))
    matrix <- matrix[h.r$order, h.c$order]
  }

  x <- x.size * 1:nrow(matrix)
  y <- y.size * 1:ncol(matrix)
  breaks <- (0:num.colors)/num.colors
  image(x, y, matrix, xlab=xlab, ylab=ylab, axes=FALSE, col=col, breaks=breaks)
  matrix
}

my.colorbar <- function(n=65, col=rg.array.colors) {
  # By default, Matlab plots n=65.
  I <- (n-1):0 / (n-1)
  M <- matrix(I, length(I), 1)
  my.heatmap(M, col=col)
}

.groupplot.calc.X <- function(
  Y, col.width=1.0, glyph.width=0.2, glyph.height=0.2, offset=0) {
  # x and y are the centers of the glyph.

  # Group each of the Y coordinates into bins of glyph_height.
  bin.height <- glyph.height
  bins <- sapply(Y, function(y) floor(y/bin.height))

  max.glyphs <- 1 + 2*col.width/glyph.width
  X.new <- c(); Y.new <- c(); I.new <- c()
  for(bin in unique(bins)) {
    I <- which(bins==bin)
    I <- I[order(Y[I])]
    Y.bin <- Y[I]

    delta <- glyph.width
    if(length(Y.bin) > max.glyphs)
      delta <- col.width / ((length(Y.bin)-1)/2)

    X.bin <- sapply(1:length(Y.bin), function(i) (i-1)/2*delta)
    if(length(X.bin) >= 2)
      X.bin[seq(2, length(X.bin), 2)] <- -X.bin[seq(2, length(X.bin), 2)]

    X.new <- c(X.new, X.bin)
    Y.new <- c(Y.new, Y.bin)
    I.new <- c(I.new, I)
  }

  X.new <- X.new + offset
  list(X.new, Y.new, I.new)
}

sample.colorbars <- function() {
  M <- t(matrix(1:8, 4, 2))
  layout(M)
  col.fns <- c(rg.array.colors, by.array.colors, rgb.colors, ryb.colors, 
    matlab.colors, matlab.hot.colors, genespring.colors, yahoo.weather.colors)
  col.names <- c("Red/Green", "Yellow/Blue", "RGB", "RYB", 
    "Matlab Jet", "Matlab Hot", "GeneSpring", "Weather")
  # For some reason, I can't loop 1:length(col.fns).  col <-
  # col.fns[i] always gets me rg.array.colors.
  i <- 1
  for(col in col.fns) {
    my.colorbar(col=col)
    title(col.names[i])
    i <- i + 1
  }
}

my.groupplot <- function(Ys, col.width=1.0, glyph.width=0.03, 
  glyph.height=0.03, col=NA) {
  # Ys is a list of vectors, where each vector is the Y-coordinates of
  # the points to plot.
  if(length(col.width) != 1) stop("invalid col.width")
  if(length(col)==1 && is.na(col))
    col <- lapply(Ys, function(y) rep("#000000", length(y)))

  num.groups <- length(Ys)

  offset = col.width
  usable.col.width = col.width/2.0 * 0.80

  Xs.new <- c(); Ys.new <- c(); col.new <- c()
  for(i in 1:num.groups) {
    Y <- Ys[[i]]
    x <- .groupplot.calc.X(
      Y, col.width=usable.col.width, glyph.width=glyph.width,
      glyph.height=glyph.height, offset=offset)
    X <- x[[1]]; Y <- x[[2]]; I <- x[[3]]
    Xs.new <- c(Xs.new, X)
    Ys.new <- c(Ys.new, Y)
    col.new <- c(col.new, col[[i]][I])
    offset <- offset + col.width
  }

  plot(Xs.new, Ys.new, pch=19, cex=1, frame.plot=FALSE, axes=FALSE, 
    col=col.new, xlab="", ylab="")
}

# Make a PCA plot of the samples (columns).
my.pcaplot <- function(X, K=3, d1=1, d2=2) {
  K <- max(c(K, d1, d2))
  S <- svd(X)
  Y <- t(S$u[,1:K]) %*% X
  plot(Y[d1,], Y[d2,], pch=".", cex=8, xlab="", ylab="")
  Y
}

my.violin <- function(DATA, 
  title="", subtitle="", xlab="", ylab="",
  text.size=20, xlab.size=20, ylab.size=20,
  xlabel.orientation="horizontal",
  log.scale=FALSE,
  show.points=FALSE, dot.size=0.5, 
  col="#000000", col.outline="#000000") {
  # DATA is a data frame with members "Category" and "Value".
  library(ggplot2)

  #lwd_box = 2
  ##lwd_axis = 2
  ##lwd_regr = 3
  #cex = 1.0 * scale_points
  #cex_lab = 1.5
  #cex_main = 2.0
  #cex_sub = 1.0
  #cex_xlab = 1.5 * scale_xlab

  p <- ggplot(DATA, aes(x=Category, y=Value))
  p <- p + geom_violin(fill=col, colour=col.outline)
  
  if(log.scale) {
    p <- p + scale_y_continuous(trans="log10", labels=sciNotation)
  }
  
  # Add dot plot.
  if(show.points) {
    p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=dot.size)
  }

  # Jittered points.
  #p <- p + geom_jitter(shape=16, position=position_jitter(0.2))
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major=element_blank())
  p <- p + theme(panel.grid.minor=element_blank())
  p <- p + theme(panel.border=element_rect(size=2))
  # Put a margin around the whole plot.
  #p <- p + theme(plot.margin=
  #  margin(mar.top, mar.right, mar.bottom, mar.left, unit="in"))
  p <- p + theme(plot.margin=
    margin(20, 20, 20, 50, unit="pt"))
  p <- p + labs(title=title, x=xlab, y=ylab) + 
    theme(plot.title=element_text(hjust=0.5))
  p <- p + labs(subtitle=subtitle)
  p <- p + theme(text=element_text(size=text.size))
  
  # Move the x-axis label away from the plot.  Negative is down.
  angle <- 0
  hjust <- 0.5
  vjust <- 1.0
  if(xlabel.orientation == "vertical") {
    angle <- 90
    hjust <- 1.0
    vjust <- 0.5
  } else if(xlabel.orientation == "diagonal") {
    angle <- 45
    hjust <- 1.0
    vjust <- 1.0
  }
  p <- p + theme(axis.text.x=element_text(
    size=xlab.size, angle=angle, hjust=hjust, vjust=vjust))
  p <- p + theme(axis.title.x=element_text(vjust=-40))
  # Move the y-axis label away from the plot.  Positive is left.
  p <- p + theme(axis.text.y=element_text(size=ylab.size))
  p <- p + theme(axis.title.y=element_text(vjust=5))
      
  # Move the x-axis tick labels away from the plot.  Doesn't work?
  #R("p <- p + theme(axis.text.x=element_text(vjust=100))")
  #R("p <- p + opts(axis.title.x=theme_text(vjust=2))")
  # Draw a line through the median.
  p <- p + stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
    geom="crossbar", size=0.5, width=0.5, color="red")
  print(p)   # stupid
}


my.beeswarm <- function(DATA, 
  # Plot
  show.box=TRUE,
  lwd.box=2,
  # Titles and labels
  title="", subtitle="", x.title="", y.title="",
  cex.title=2.0,
  cex.subtitle=1.0,
  cex.axis.title=1.5,
  cex.x.label=1.5,
  cex.y.label=1.5,
  xlabel.orientation="horizontal",
  ylabel.orientation="horizontal",
  # Axes
  ymin=NULL,
  ymax=NULL,
  # Points
  cex.points=1.0,
  color="#000000",
  log.scale=FALSE,
  plot.means=FALSE) {
  # DATA is a list of vectors.  The names are the names of the
  # categories, and the vectors contain the data points.
  library(beeswarm)

  ylim <- NULL
  if(!is.null(ymin) || !is.null(ymax)) {
    x <- c()
    for(i in 1:length(DATA))
      x <- c(x, DATA[[i]])
    if(is.null(ymin))
      ymin <- min(x)
    if(is.null(ymax))
      ymax <- max(x)
    if(log.scale)
      ymin <- max(ymin, 1)
    ylim <- c(ymin, ymax)
  }

  if(log.scale) {
    # Remove non-positive values.
    for(i in 1:length(DATA))
      DATA[[i]] <- DATA[[i]][DATA[[i]]>0]
  }

  #lwd_axis = 2
  #lwd_regr = 3
  #cex_names = 1.5 * scale_labels

  # cex.axis.title controls the size of both the X- and Y- axis titles

  x <- beeswarm(DATA, main="", xlab="", ylab="", pch=19, cex=cex.points,
    ylim=ylim, corral="wrap", spacing=1.0, log=log.scale, axes=FALSE,
    col=color)

  if(plot.means) {
    means <- unlink(lapply(DATA, mean))
    for(i in 1:length(means)) {
      segments(i+0.8, means[i], i+1.2, means[i], lwd=2, lty="solid", col="red")
    }    
  }

  # Make plot area solid white.
  #jmath.R('usr <- par("usr")')
  #jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
  #jmath.R_fn(
  #    "hist", jmath.R_var("X"), plot=jmath.R_var("FALSE"),
  #    main=main, xlab="", ylab="", axes=jmath.R_var("FALSE"),
  #    add=jmath.R_var("TRUE"))

  # Calculate correlation, and other statistics.
  # TODO: Should calculate this for each series.
  #r = jmath.R("cor(X, Y)")
  #p_value = jmath.R("cor.test(X, Y)$p.value")
  #r = r[0]
  #p_value = p_value[0]
  #print "R = %.2f" % r
  #print "p = %.2g" % p_value

  if(show.box)
    box(lwd=lwd.box)

  # X-axis
  x.labels <- names(DATA)
  x.at <- 1:length(x.labels)
  
  if(xlabel.orientation == "vertical") {
    las <- 3   # vertical labels
    axis(1, labels=x.labels, at=x.at, las=las, cex.axis=cex.x.label)
  } else if(xlabel.orientation == "diagonal") {
    y.line.height <- par("cxy")[2]
    axis(1, labels=FALSE, at=x.at)
    y <- par("usr")[3]    # bottom of plot, in user coordinates
    # Given in logged scale.  Have to unlog for text.
    if(log.scale) {
      y <- 10 ** y
      # Go down 1 line in user coordinates.
      y <- y - 0.10
    } else {
      # Go down 1 line in user coordinates.
      y <- y - y.line.height * 1
    }
    # When log.scale is TRUE, then this is in logged coordinates.
    # xpd=TRUE allows drawing of labels outside clipping region.
    text(x=x.at, y=y, labels=x.labels, srt=45, adj=1, 
      cex=cex.x.label, xpd=TRUE)
    # xps=TRUE, 
  } else {
    #las <- 1   # horizontal labels
    axis(1, labels=x.labels, at=x.at, cex.axis=cex.x.label)
  }

  # Y-axis
  y.ticks <- axTicks(2, log=log.scale)
  y.at <- y.ticks
  y.labels <- y.ticks
  if(log.scale) {
    exp.min <- min(floor(log(y.ticks, 10)))
    exp.max <- max(floor(log(y.ticks, 10)))
    y.at <- outer(1:9, 10^(exp.min:exp.max))
    y.lab <- ifelse(log10(y.at) %% 1 == 0, y.at, NA)
    y.lab <- y.lab[y.at <= max(y.ticks)]
    y.at <- y.at[y.at <= max(y.ticks)]
    x <- sciNotation(y.lab, 1)
    x[is.na(y.lab)] <- NA
    y.labels <- x
  }

  y.las <- 0  # parallel to axis, default
  if(ylabel.orientation == "vertical")
    y.las <- 1
  axis(2, labels=y.labels, at=y.at, las=y.las, cex.axis=cex.y.label)

  title(main=title, sub=subtitle, xlab=x.title, ylab=y.title,
    cex.lab=cex.axis.title, cex.main=cex.title, cex.sub=cex.subtitle)
}

my.boxplot <- function(X, labels=TRUE, col=NULL, main="", xlab="", ylab="", 
  leg=NULL, fill=NULL, cex.labels=1.25, cex.legend=1, vert.labels=TRUE,
  ylim=NULL, cex.lab=1.5, sub="", cex.sub=1.5, x.legend=NULL)
{
  las <- 0
  if(vert.labels)
    las <- 3
  at <- NULL
  if(length(labels) > 1)
    at <- 1:length(labels)
  if(is.null(x.legend))
    x.legend <- "bottomleft"

  lwd <- 2
  boxplot(X, col=col, xlab="", ylab="", axes=FALSE, pch=19, cex=1, ylim=ylim)
  box(lwd=lwd)
  axis(1, lwd=lwd, cex.axis=cex.labels, labels=labels, at=at, las=las)
  axis(2, lwd=lwd, cex.axis=1.5)
  title(main=main, xlab=xlab, ylab=ylab, cex.lab=cex.lab, cex.main=2.0, 
    sub=sub, cex.sub=cex.sub)
  if(length(leg))
    legend(x.legend, legend=leg, fill=fill, box.lwd=1.5, cex=cex.legend,
      bg="#FFFFFF", inset=0.05)
}


my.barplot <- function(DATA, bar.names, 
    # Titles and labels.
    title="", subtitle="", x.title="", y.title="",
    cex.title=2.0, cex.subtitle=1.0, 
    cex.x.title=1.5, cex.y.title=2,
    # Axes
    ymin=NULL, ymax=NULL,
    lwd.axis=2,
    show.x.axis=TRUE, show.y.axis=TRUE,
    show.box=TRUE,
    # Legend
    show.legend=FALSE,
    legend.names=NULL,
    legend.loc="topright",
    legend.inset=0.05,
    lwd.legend=2,
    cex.legend=1.5,  # How big to make the legend
    color.legend=NULL,
    # Bars
    label.bars=TRUE,
    cex.bar.lab=1.5,
    bar.label.orientation="horizontal",
    show.frequencies=FALSE,
    density=NULL,
    lwd.box=4, color="#000000"
    ) {

    ylim <- NULL
    if(!is.null(ymin) & !is.null(ymax)) {
        ylim <- c(ymin, ymax)
    } else if(!is.null(ymin)) {
        ylim <- c(ymin, max(DATA))
    } else if(!is.null(ymax)) {
        ylim <- c(min(DATA), ymax)
    }

    if(show.frequencies) {
        # Add space for two lines of text.
        y.per.line <- par("cxy")[2]
        ylim[2] <- ylim[2] + y.per.line*2
    }

    beside <- FALSE
    if(nrow(DATA) > 1)
        beside <- TRUE

    # Need to do this, or else colors won't work.
    if(nrow(DATA) == 1)
        DATA <- as.numeric(DATA)
    x <- barplot(DATA, beside=beside, main="", xlab="", ylab="",
        xpd=FALSE,   # so bars don't bleed
        col=color,
        density=density,
        axes=FALSE,
        ylim=ylim,
        yaxt="n",   # don't print y axis
        xaxt="n",   # don't print x axis
        names.arg=NULL
        )

    # Make vectors, in the order of the bars.
    X.coords <- as.numeric(x)
    Y.coords <- as.numeric(values)
    
    if(show.frequencies) {
        text(x=X.coords, y=Y.coords, label=Y.coords, pos=3, cex=1.0, col="red")
    }

    if(show.legend & (length(color) > 1)) {
        if(length(legend.names) != length(color)) stop("legend names")
        legend(legend.loc, legend=legend.names, fill=color,
               inset=legend.inset, cex=cex.legend, bty="n",
               box.lwd=lwd.legend)
    }
    if(!is.null(color.legend)) {
        # Parse color_legend.  Should be list of "#FFFFFF,<name>"
        colors <- c()
        names <- c()
        for(i in 1:length(color.legend)) {
            y <- strsplit(color.legend[i], ",")[[1]]
            c <- y[1]
            n <- paste(y[2:length(y)], collapse=",")
            colors <- c(colors, c)
            names <- c(names, n)
        }
        legend(legend.loc, legend=names, fill=colors,
            inset=legend.inset, cex=cex.legend, bty="n", box.lwd=lwd.legend)
    }

    if(show.box)
        box(bty="l", lwd=lwd.box)

    par(font.lab=2)
    
    at <- NULL

    if(label.bars) {
        if(nrow(values) > 1) {
            # X.coords is in order:
            #   series0 value0, series1 value0, series2 value0,
            #   series0 value1, series1 value1, series2 value1, etc.
            # at should be in order:
            #   value0, value1, value2
            num.series <- nrow(values)
            at = c()
            for(i in seq(1, length(X.coords), num.series)) {
                x <- X.coords[i:(i+num.series-1)]
                at <- c(at, mean(x))
                }
        } else {
            at <- X.coords
        }
    }

    if(!show.x.axis | !label.bars) {
    } else if(bar.label.orientation == "vertical") {
        las <- 3   # vertical labels
        axis(1, lwd=lwd.axis, xpd=NA, at=at, labels=bar.names,
           tick=FALSE, las=las, cex.axis=cex.bar.lab)
    } else if(bar.label.orientation == "diagonal") {
        y.line.height <- par("cxy")[2]
        x <- at
        y <- par("usr")[3]    # bottom of plot, in user coordinates
        # Go down 1 line in user coordinates.
        y <- y - y.line.height * 1
        text(x=x, y=y, labels=bar.names, srt=45, adj=1, xpd=NA, 
            cex=cex.bar.lab)
    } else {
        las <- 1   # horizontal labels
        axis(1, lwd=lwd.axis, xpd=NA, at=at, labels=bar.names, tick=FALSE,
            las=las, cex.axis=cex.bar.lab)
    }
    if(show.x.axis) {
        mtext(x.title, side=1, line=2.8, cex=cex.x.title)
    }
    if(show.y.axis) {
        axis(2, lwd=lwd.axis, xpd=NA, cex.axis=cex.y.title)
        mtext(y.title, side=2, line=2.8, cex=cex.y.title)
    }
  
    title(main=title, sub=subtitle, cex.main=cex.title, cex.sub=cex.subtitle)
}


# Given the name of a file, try to figure out the type of bitmap to
# create.
get.bitmap.type <- function(filename) {
  # Default is PNG.
  bm.type <- "png16m"
  x <- tolower(filename)
  x <- substring(x, nchar(x)-3, nchar(x))
  if(x == ".pdf")
    bm.type <- "pdfwrite"
  bm.type
}

sciNotation <- function(x, digits=1) {
  if(length(x) > 1)
    return(append(sciNotation(x[1]), sciNotation(x[-1])))
  if(is.na(x) | !x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  if(base == 1)
    return(as.expression(substitute(10^exponent,
      list(exponent=exponent))))
  as.expression(substitute(base %*% 10^exponent,
    list(base=base, exponent=exponent)))
}

calc.plot.margins <- function(
  temp.file, desired.width, desired.height, 
  pointsize, x.labels=NULL, x.orientation="horizontal", x.cex=1,
  x.title.cex=NULL, y.labels=NULL, y.cex=1, y.title.cex=NULL) {
  # Need to create a plot so that we can get the default plotting
  # area, and calculate the size of the text.
  # temp.file       temporary plot file
  # desired.width   width, in inches
  # desired.height  height, in inches
  # pointsize       point size for rendering text
  # x.labels        vector of x labels
  # x.orientation   "horizontal", "vertical", "diagonal".  partial match
  # x.cex           size of X-axis labels
  # x.title.cex     size of X-axis title.  NULL if no title.
  # y.labels
  # y.title.cex     size of Y-axis title.  NULL if no title.
  # 
  # Returns a list with members:
  # height        total height to make plot (all are in inches)
  # width         total width to make plot
  # pointsize
  # plot.height   height of plotting area (not including margins)
  # plot.width    width of plotting area (no margins)
  # mai           margins, in inches
  # x.label.line
  # x.title.line
  # y.label.line
  # y.title.line
  #
  # This works for standard R.  It will only provide an estimate for ggplot.

  ORIENTATION <- c("vertical", "horizontal", "diagonal")
  x.orient <- pmatch(x.orientation, ORIENTATION)
  if(is.na(x.orient)) stop("unknown orientation")

  bitmap(temp.file, type="pdfwrite", res=300, pointsize=pointsize,
    height=desired.height, width=desired.width, units="in")
  plot.new()

  # Calculate the size of the plotting area (not including margins).
  x <- par("pin")
  plot.width <- x[1]
  plot.height <- x[2]


  # Calculate the size of the bottom margin.
  line.height <- par("csi")    # will create devices.  don't close yet

  X.LAB.LINE <- 1
  x.label.line <- X.LAB.LINE
  x.label.width <- 0
  x.label.height <- 0
  if(!is.null(x.labels)) {
    # Calculate the sizes of the x.labels.
    widths <- sapply(x.labels, function(x) 
      strwidth(x, units="inches", cex=x.cex))
    heights <- sapply(x.labels, function(x) 
      strheight(x, units="inches", cex=x.cex))

    # If the labels are vertical, then rotate them.
    if(x.orient == 1) {          # vertical
      x <- widths
      widths <- heights
      heights <- x
    } else if(x.orient == 3) {   # diagonal
      # Calculate for a 45 degree counter-clockwise rotation.
      # Assume no height.  May be off a bit.
      x1 <- sin(45) * widths
      x2 <- cos(45) * widths
      heights <- x1
      widths <- x2
    }
    x.label.width <- max(widths)
    x.label.height <- max(heights)
  }

  x.title.height <- 0
  if(!is.null(x.title.cex))
    x.title.height <- strheight("W", units="inches", cex=x.title.cex)
  # Add a boundary between labels and title.
  x.title.line <- x.label.line + (x.label.height*1.5/line.height)

  # Calculate the number of inches available for the X-axis tick
  # labels.  This is useful if the x-axis needs to be vertical.  Since
  # the X-axis tick labels are printed on line 1, inches available is
  # the total size of the bottom margin, minus the size of one line.

  # Calculate the amount of margin needed on the bottom.
  x1 <- line.height*X.LAB.LINE
  x2 <- x.label.line*line.height + x.label.height
  x3 <- x.title.line*line.height + x.title.height
  # Add a small buffer.
  bottom.needed <- max(c(x1, x2, x3)) + line.height

  # Calculate how much to increase the bottom margin.
  mar.bottom <- par("mai")[1]
  add.to.bottom <- max(bottom.needed - mar.bottom, 0)


  # Calculate the size of the left margin.
  Y.LAB.LINE <- 1
  y.label.line <- Y.LAB.LINE
  y.label.height <- 0
  if(!is.null(y.labels)) {
    heights <- sapply(y.labels, function(x) 
      strheight(x, units="inches", cex=y.cex))
    y.label.height <- max(heights)
  }
  y.title.height <- 0
  if(!is.null(y.title.cex))
    y.title.height <- strheight("W", units="inches", cex=y.title.cex)
  # Add a boundary between labels and title.
  y.title.line <- y.label.line + (y.label.height*1.5/line.height)


  # Calculate the amount of margin needed on the left.
  x1 <- line.height*Y.LAB.LINE
  x2 <- y.label.line*line.height + y.label.height
  x3 <- y.title.line*line.height + y.title.height
  # Add a small buffer.
  left.needed <- max(c(x1, x2, x3)) + line.height

  # Calculate how much to increase the left margin.
  mar.left <- par("mai")[2]
  add.to.left <- max(left.needed - mar.left, 0)

  # Calculate the final extents for the figure.
  height <- desired.height + add.to.bottom
  width <- desired.width + add.to.left

  mai.bottom <- par("mai")[1] + add.to.bottom
  mai.left <- par("mai")[2] + add.to.left
  mai.top <- par("mai")[3]
  mai.right <- par("mai")[4]
  mai <- c(mai.bottom, mai.left, mai.top, mai.right)
  
  # Don't need plot anymore.  Turn off.
  dev.off()
  unlink(temp.file)

  x <- list(height=height, width=width, pointsize=pointsize,
    plot.height=plot.height, plot.width=plot.width, mai=mai,
    x.label.line=x.label.line, x.title.line=x.title.line, 
    y.label.line=y.label.line, y.title.line=y.title.line
    )
  x
}

normalize.one.mv <- function(x, M=0, V=1) {
  # Normalize a list of numbers so the mean is M and variance is V.
  M.0 <- mean(x)
  V.0 <- var(x)
  if(is.null(M))
    M <- M.0
  if(is.null(V))
    V <- V.0
  if(V.0 == 0)
    return(x-M.0+M)
  (x-M.0)*sqrt(V/V.0) + M
}

normalize.mv <- function(X, M=0, V=1) {
  t(apply(X, 1, function(x) normalize.one.mv(x, M, V)))
}

matrix2color <- function(cmatrix, pos) {
  # pos is [0, 1].  Returns r, g, b where each one is from [0, 1].
  if(is.nan(pos))  # this can happen if someone calls matlab.colors(1)
    pos <- 0.5
  breaks <- cmatrix[,1]
  i1 <- sum(pos >= breaks)
  x <- cmatrix[i1,2:4]
  if(i1 < nrow(cmatrix)) {
    i2 <- i1 + 1
    delta <- (pos - cmatrix[i1,1]) / (cmatrix[i2,1]-cmatrix[i1,1])
    x <- cmatrix[i1,2:4] + delta*(cmatrix[i2,2:4]-cmatrix[i1,2:4])
  }
  x/255
}

matrix2rgb <- function(cmatrix, pos) {
  # pos is [0, 1].  Returns color, e.g. "#003300".
  x <- matrix2color(cmatrix, pos)
  rgb(x[1], x[2], x[3], maxColorValue=1)
}

