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

