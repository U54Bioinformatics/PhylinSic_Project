library(argparse)

parser <- ArgumentParser()

parser$add_argument("--input_file", "-i", 
  help="File of sequences to model in FASTA format.")
parser$add_argument("--output_dir", "-o", help="Output directory")

parser$add_argument("--beast2_dir")
parser$add_argument("--site_model")
parser$add_argument("--clock_model")
parser$add_argument("--tree_prior")
parser$add_argument("--iterations", type="integer")
parser$add_argument("--sample_interval", type="integer")
parser$add_argument("--rng_seed", type="integer")

args <- parser$parse_args()


#FASTA.FILE <- snakemake@input[[1]]
#OUT.DIR <- snakemake@output[[1]]
FASTA.FILE <- args$input_file
OUT.DIR <- args$output_dir

BEAST2.DIR <- args$beast2_dir
BEAST.INFILE <- file.path(OUT.DIR, "beast2.infile.txt")
BEAST.OUTFILE <- file.path(OUT.DIR, "beast2.outfile.txt")
BEAST.MODEL.FILE <- file.path(OUT.DIR, "beast2.model.RDS")
TREE.LOG <- file.path(OUT.DIR, "tree.log")
TRACE.LOG <- file.path(OUT.DIR, "trace.log")
SCREEN.LOG <- file.path(OUT.DIR, "screen.log")

#SITE.MODEL <- snakemake@params$site_model
#CLOCK.MODEL <- snakemake@params$clock_model
#TREE.PRIOR <- snakemake@params$tree_prior
#NUM.ITER <- snakemake@params$iterations
#SAMPLE.INTERVAL <- snakemake@params$sample_interval
#RNG.SEED <- snakemake@params$rng_seed

SITE.MODEL <- args$site_model
CLOCK.MODEL <- args$clock_model
TREE.PRIOR <- args$tree_prior
NUM.ITER <- args$iterations
SAMPLE.INTERVAL <- args$sample_interval
RNG.SEED <- args$rng_seed

#source("/usr/local/changlab/Rlib/filelib.R")
#source("/usr/local/changlab/Rlib/plotlib.R")
source("scripts/beastlib.R")
library(parallel)

library(ggplot2)
library(babette)
library(ape)
library(phangorn)
library(TreeTools)

library(dplyr)
library(data.table)
library(ggtree)

#if(!dir.exists(OUT.DIR))
#  dir.create(OUT.DIR, recursive=TRUE, mode="0755")

# Read the FASTA file.
x <- read.phyDat(FASTA.FILE, format="fasta")
# Names of cells cannot contain punctuation other than "_" ... this
# is a check and substitution for that.
names(x) <- gsub("[^[:alnum:]\\\\-\\\\s]", "_", names(x))
# cell x position matrix.  rownames are the cell names.
M.seqs <- PhyDatToMatrix(x)

print(sprintf("Start with %d sequences and %d positions.", nrow(M.seqs),
  ncol(M.seqs)))
if(min(dim(M.seqs)) < 2) stop("Not enough data.")


sparsity <- sum(M.seqs == "-")/length(M.seqs)*100
print(sprintf("Final data has %d sequences and %d positions, %.1f%% sparse.",
  nrow(M.seqs), ncol(M.seqs), sparsity))


## Save filtered data in fasta format.
#x <- MatrixToPhyDat(M.seqs)
#write.phyDat(x, file=FILTERED.FASTA.FILE, format="fasta")

# Fit the data.
x <- formatC(NUM.ITER, format="f", big.mark=",", digits=0)
print(sprintf("Fitting tree over %s iterations.", x))




# Run BEAST.
launcher.path <- sprintf("%s/lib/launcher.jar", BEAST2.DIR)
if(!file.exists(launcher.path)) stop("no launcher.jar")

start.time <- Sys.time()
beast.model <- run.beast(FASTA.FILE, SITE.MODEL=SITE.MODEL,
  CLOCK.MODEL=CLOCK.MODEL, TREE.PRIOR=TREE.PRIOR, num.iter=NUM.ITER,
  sample.interval=SAMPLE.INTERVAL, rng.seed=RNG.SEED,
  beast2.path=launcher.path,
  beast.infile=BEAST.INFILE, beast.outfile=BEAST.OUTFILE,
  tree.log=TREE.LOG, trace.log=TRACE.LOG, screen.log=SCREEN.LOG)
end.time <- Sys.time()
secs <- as.numeric(end.time) - as.numeric(start.time)
print(sprintf("run time %.3fs", secs))



total.iters <- max(beast.model$estimates$Sample)
print(sprintf("total %g iterations", total.iters))


# Write out the BBT model.
saveRDS(beast.model, file=BEAST.MODEL.FILE)

print("The beast is done")
