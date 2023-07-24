library(argparse)

parser <- ArgumentParser()
parser$add_argument("filenames", nargs="+")
args <- parser$parse_args()

MODEL.FILE <- args$filenames[1]
I.MIN <- as.numeric(args$filenames[2])
I.MAX <- as.numeric(args$filenames[3])
POSTERIOR.FILE <- args$filenames[4]
POSTERIOR.ZOOM.FILE <- args$filenames[5]
TREE.HEIGHT.FILE <- args$filenames[6]
YULE.MODEL.FILE <- args$filenames[7]
NBURNIN <- args$filenames[8]

library(ggplot2)

XLAB.SIZE <- 32
YLAB.SIZE <- 32
XTICK.SIZE <- 24
YTICK.SIZE <- 24

VLINE.SIZE <- 2
LINE.SIZE <- 2

my.theme <- theme(
  axis.title.x=element_text(size=XLAB.SIZE, margin=margin(t=12,r=0,b=0,l=0)),
  axis.title.y=element_text(size=YLAB.SIZE, margin=margin(t=0,r=12,b=0,l=0)),
  axis.text.x=element_text(size=XTICK.SIZE),
  axis.text.y=element_text(size=YTICK.SIZE)
  )

beast.model <- readRDS(MODEL.FILE)

if(is.null(beast.model$estimates$Sample))
  stop("missing Sample from beast model file")

# Plot the posterior.
p <- ggplot(data=beast.model$estimates) +
  theme_classic() + my.theme + 
  xlab("MCMC Iteration") +
  ylab("Posterior")
if(!is.null(NBURNIN))
  p <- p + geom_vline(xintercept=NBURNIN, linetype="dotted", size=VLINE.SIZE)
p <- p + geom_line(
  aes_string(x="Sample", y="posterior"), color="red", size=LINE.SIZE)

pdf(POSTERIOR.FILE, height=8, width=11)
print(p)
dev.off()


# Plot the zoomed in posterior.
if(I.MIN < 1) stop("bad 1")
if(I.MAX > nrow(beast.model$estimates)) stop(
  sprintf("bad 2 %s %s", I.MAX, nrow(beast.model$estimates)))
if(I.MIN >= I.MAX) stop("bad 3")
p <- ggplot(data=beast.model$estimates[I.MIN:I.MAX,]) +
  theme_classic() + my.theme + 
  xlab("MCMC Iteration") +
  ylab("Posterior")
p <- p + geom_line(
  aes_string(x="Sample", y="posterior"), color="red", size=LINE.SIZE)

pdf(POSTERIOR.ZOOM.FILE, height=8, width=11)
print(p)
dev.off()



# Plot the tree height.
p <- ggplot(data=beast.model$estimates) +
  theme_classic() + my.theme + 
  xlab("MCMC Iteration") +
  ylab("Tree Height")
if(!is.null(NBURNIN))
  p <- p + geom_vline(xintercept=NBURNIN, linetype="dotted", size=VLINE.SIZE)
p <- p + geom_line(
  aes_string(x="Sample", y="TreeHeight"), color="green", size=LINE.SIZE)

pdf(TREE.HEIGHT.FILE, height=8, width=11)
print(p)
dev.off()


# Plot the YuleModel, if given.
if(!is.null(beast.model$estimates$YuleModel)) {
  p <- ggplot(data=beast.model$estimates) +
    theme_classic() + my.theme +
    xlab("MCMC Iteration") +
    ylab("Yule Model")
  if(!is.null(NBURNIN))
    p <- p + geom_vline(xintercept=NBURNIN, linetype="dotted", size=VLINE.SIZE)
  p <- p + geom_line(
    aes_string(x="Sample", y="YuleModel"), color="blue", size=LINE.SIZE)

  pdf(YULE.MODEL.FILE, height=8, width=11)
  print(p)
  dev.off()
}
