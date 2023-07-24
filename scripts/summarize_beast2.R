library(argparse)

parser <- ArgumentParser()
parser$add_argument("filenames", nargs="+")
args <- parser$parse_args()

MODEL.FILE <- args$filenames[1]
SAMPLE.INTERVAL <- as.numeric(args$filenames[2])
NBURNIN <- as.numeric(args$filenames[3])
OUTFILE <- args$filenames[4]
ESS.FILE <- args$filenames[5]



library(tracerer)
source("scripts/filelib.R")

beast.model <- readRDS(MODEL.FILE)

total <- max(beast.model$estimates$Sample)
if(total <= 0) stop("no samples")
if(NBURNIN >= total) stop("bad: burn in >= iterations")
if(NBURNIN >= (total-SAMPLE.INTERVAL))
  stop("Too many iterations being discarded from burn-in")
frac <- NBURNIN / total

traces <- remove_burn_ins(traces=beast.model$estimates,
  burn_in_fraction=frac)

# 2 is OK.  No problems.  Not sure about 1.
if(length(traces$posterior) < 2)
  stop("Burn in too high.  No more iterations after removing burn in.")

# Write summary table.
stats <- t(
  calc_summary_stats(traces$posterior, sample_interval=SAMPLE.INTERVAL))
colnames(stats) <- "Statistic"
data.out <- cbind(rownames(stats), stats)
colnames(data.out)[1] <- "Name"
my.write(data.out, OUTFILE, col.names=TRUE)

# Write ESS table.
esses <- t(calc_esses(traces, sample_interval=SAMPLE.INTERVAL))
colnames(esses) <- "ESS"
data.out <- cbind(rownames(esses), esses)
colnames(data.out)[1] <- "Name"
my.write(data.out, ESS.FILE, col.names=TRUE)
