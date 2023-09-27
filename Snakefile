# This is a Snakefile that contains rules that implement
# phylogenetic inference from single-cell RNA-Seq data as
# described in:
#   Liu X, Griffiths J, Bishara I, Liu J, Bild AH, and Chang JT.
#   Phylogenetic inference from single-cell RNA-seq data.
#   bioRxiv 2022.09.27.509725.  doi:
#   https://doi.org/10.1101/2022.09.27.509725
#
# It should be executed with Snakemake.
#
# There are parameters below that you can customize for your
# analysis.

# Versions:
#  1  230721  Initial release.
#  2  230921  Miscellaneous bug fixes, multiple batches.
VERSION = 2



# PARAMETERS - CRITICAL
# ---------------------


# This pipeline requires the following software.  Set these parameters
# to the commands used to execute them.

SAMTOOLS = "samtools"
PICARD = "picard"
GATK = "gatk"
VARSCAN = "varscan"

# In the last step of the pipeline, we run phylogenetic inference
# by calling BEAST2 (a Java program) from R using the babette
# library.  While all the software required in the rest of the
# pipeline could be installed correctly in a Conda environment, I
# had difficulties with this step (see INSTALLATION section in
# the README.md file).  Therefore, we installed this
# software (Java, R+babette, BEAST2) outside Conda.  These
# parameters can be used to customize how this software is run.
#
# If BEAST2+dependencies can run from Conda correctly on your
# system (congratulations!), you do not need to change these
# parameters.  Otherwise, point them to your working
# installation.
#
# If you're not sure what to do, install this software from Conda
# first.  If it doesn't work, then either try to fix the Conda
# installation, or install it yourself separately.

RSCRIPT = "Rscript"
JAVA = "java"

# This should point to the directory that beast is installed in, e.g.
#   /usr/local/beast
#
# If set to None, I will try to find it myself.
BEAST2_DIR = None

# This should point to the "logcombiner" program that is installed with
# BEAST2.
LOGCOMBINER = "logcombiner"






# PARAMETERS - VARIANT CALLING
# ----------------------------
# We call variants on a "pseudobulk" sample where the reads from
# all cells are merged.  These parameters are used to select for
# high quality variants.


# Each of these can filters can be disabled by setting to None.

# Only keep the variants support by at least this many alternate
# reads.
KEEP_CALLS_WITH_MIN_ALT_READS = 5

# Only keep the variants support by at least this many total
# reads.
KEEP_CALLS_WITH_MIN_TOTAL_READS = 20

# Only keep the variants support by at least this variant allele
# frequency.
KEEP_CALLS_WITH_MIN_VAF = 0.05



# PARAMETERS - SELECTING SITES
# ----------------------------
# The phylogeny is generated on only a subset of the variant
# sites because 1) it takes too much computation to use all the
# potential sites, and 2) most of the variants from single-cell
# RNA-Seq are poor anyway.  We select high quality sites using a
# range of criteria.

# Each of these can filters can be disabled by setting to None.

# Only keep the sites if at least 1 cell has this number of reads.
KEEP_SITES_WITH_HIGH_READS = 5

# Remove the sites that reside on chromosomes that start with one
# of these prefixes.  This should be a comma-separated list.
REMOVE_CHROM_WITH_PREFIX = "MT,Y,GL"

# Remove the sites that have more than one alternate allele,
# e.g. A -> C or G.  This should be "yes" or "no".
REMOVE_MULTIALLELIC_SITES = "yes"

# Remove sites that are < this number of bases apart.  This can
# be indicative of an alignment error or high discrepancy region.
# If two sites are very close together, arbitrarily drop one of
# them.
REMOVE_CLOSE_SITES = 5

# Select sites that exhibit variation in the genotype across
# cells.  If 100% of the sites have the same genotype across all
# cells, then it's not informative for generating a
# phylogeny (because there is no difference in the genotype at
# this site).  We can select for the sites that have the most
# variation by keeping the ones where the minor allele genotype
# is seen in >= X% of the cells.  This parameter controls X.
KEEP_SITES_WITH_MIXED_GENOTYPES = 20


# Select the top N sites, where the "top" site is the one seen in
# the most cells (i.e. not dropout).
KEEP_N_SITES_SEEN_IN_MOST_CELLS = 1000



# PARAMETERS - GENOTYPE CALLS AND IMPUTATION
# ------------------------------------------
# These parameters control how the genotype calls are made.
# Imputation can be enabled or disabled here.

# Number of sites to use when comparing cells in kNN comparison.

# We use a nearest neighbors algorithm to smooth genotype calls.
# For each cell, we calculate its neighbors by comparing the
# distribution of reads across N sites with high coverage.  This
# controls N.
TERNARY_FEATURES = 200

# This is the number of neighbors to use to smooth the genotype
# calls.
TERNARY_K = 10

# We compare two cells by comparing sequences drawn from
# distributions determined by their read counts and scoring the
# similarity of their sequences.  This is the number of samples
# to generate.  Higher numbers will lead to more accurate scores,
# but takes longer.
TERNARY_SAMPLES = 100

# We smooth the genotypes for each cell by taking a weighted
# average of the probability distribution of the genotypes of a
# cell and its neighbors.  This controls how much weight to give
# to its neighbors.  The default value will give each cell the
# same weight.  Set to 0 to turn off smoothing.
TERNARY_DELTA = float(TERNARY_K)/(TERNARY_K+1)

# If a cells has no reads at a site, this controls whether to
# impute the genotype at that site, or to leave it as a missing
# value.
TERNARY_IMPUTE = True

# Set the seed of the random number generator.  Can use to
# make sure results are reproducible.
PHYLO_RNG_SEED = 1




# PARAMETERS - PHYLOGENETIC INFERENCE
# -----------------------------------

# Which model to use for nucleotide substitutions.  Must be one of:
# jc69    equal mutation rate
# hky     has transitions and transversions
# tn93    has different transitions
# gtr     most general
BEAST2_SITE_MODEL = "gtr"

# How to model rates of mutations across clades.  Must be one of:
# strict  same mutation rate
# rln     relaxed
BEAST2_CLOCK_MODEL = "rln"

# How to model branching in the tree.  Must be one of:
# bd      birth-death
# ccp     Coalescent Constant Population
# cep     Coalescent Exponential Population
# cbd     Coalescent Bayesian Skyline
# yule    constant birth rate
BEAST2_TREE_PRIOR = "yule"

# How many iterations to run the analysis.  For our data, we run
# ~100 million iterations or so for the final analysis, and
# shorter (e.g. 1 million) when testing.
#BEAST2_ITERATIONS = 100000000
BEAST2_ITERATIONS = 100000


# How many iterations to discard for burn-in.
BEAST2_BURNIN = int(BEAST2_ITERATIONS*0.50)


assert BEAST2_BURNIN < BEAST2_ITERATIONS, \
    "Cannot discard all samples as burnin (%d)." % BEAST2_BURNIN



# How frequently (in number of iterations) to collect statistics
# on the sampling.
BEAST2_SAMPLE_INTERVAL = 1000

# Set the random number generator seed for the tree building.
BEAST2_RNG_SEED = 1




# PARAMETERS - DATA PROCESSING
# ----------------------------
# In this pipeline, there are several steps where we run small
# batches of the data in parallel (i.e. a scatter-gather
# approach).  This can greatly decrease computation time on
# systems with multiple processors.  These parameters control the
# sizes of the batches of data.  These batch sizes are set
# relative to the number of cells that you start with.
#
# The defaults here should work for most situations, unless you
# are running on a computer that is really short on RAM.  If you
# are running into memory errors in a step that is run in
# parallel, you may want to configure a smaller batch size for
# that step.


# We demultiplex the alignments from CellRanger into smaller
# batches of cells that are preprocessed in parallel.  This
# parameter controls the number of cells in each batch.
#
# The main limitation is this parameter is the number of file
# handles that the OS can open simultaneously during the
# demultiplexing step.  Setting this number too high can exceed
# this limit, and you may get error messages that there are too
# many file handles open.
CELLS_PER_BATCH_PREPROCESSING = 16


# Calling variants from a large pseudobulk file can be very slow.
# To speed up, we split up the genome into smaller intervals and
# call variants on these intervals in parallel.  This controls
# the number of intervals used.
CELLS_PER_GENOME_INTERVAL = 64


# After calling the pseudobulk sample, we merge all the mutations
# into a single variant matrix.  This matrix may be too big to
# fit into memory, so we do it in smaller batches.
CELLS_PER_BATCH_VARIANT_MATRIX = 32


# The variant calling will almost certainly identify too many
# sites that can be processed by the phylogenetic inference
# algorithm.  Therefore, we select the highest quality sites to
# analyze and filter out the poor quality sites.  To make the
# filter run quicky and uess less memory, we filter the cells in
# parallel when possible.  This parameter controls the number of
# batches to use to filter the cells.
CELLS_PER_BATCH_FILTER1 = 512


# We smooth genotypes using a nearest neighbors algorithm.  To
# speed up the analysis, we identify the neighbors of batches of
# cells in parallel.  This controls the size of those batches.
CELLS_PER_BATCH_KNN = 512


# We call genotypes in batches of cells.  This controls the size
# of those batches.
CELLS_PER_BATCH_GENOTYPE_CALLING = 512



############################################################
# Directory Structure

# output/
#   # 1.  Prepare the reference files.
#   genome.fa
#   genome.fa.fai
#   genome.dict
#   # 2.  Demultiplex CellRanger output into single cells.
#   cells.txt             # A cleaned up version of the cells file.
#   02_demux/{sample}.bam
#   02_demux/{sample}.barcodes.txt
#   02_demux/{sample}.header.txt
#   02_demux/{sample}.alignments.txt
#   02_demux/{sample}.{batch}.barcodes.txt
#   02_demux/{sample}.{batch}.alignments.txt
#   02_demux/{sample}.{batch}.cells.sam/
#   # 3.  Preprocess the single cell files.
#   03_preproc/{sample}.{batch}.cells.bam/
#   03_preproc/{sample}.{batch}.cells.rg.bam/
#   03_preproc/{sample}.{batch}.cells.rg_contig.bam/
#   03_preproc/{sample}.{batch}.cells.rg_contig_index.bam/
#   03_preproc/{sample}.{batch}.cells.rg_contig_index_snt.bam/
#   03_preproc/{sample}.{batch}.cells.rg_contig_index_snt.report/
#   03_preproc/{sample}.{batch}.cells.rg_contig_index_snt_recal.bam/
#   # 4.  Merge single cells to a pseudobulk sample.
#   04_pbulk/{sample}.{batch}.merged.bam
#   04_pbulk/{sample}.{batch}.bam
#   04_pbulk/{sample}.merged.bam
#   04_pbulk/{sample}.bam
#   04_pbulk/pseudobulk.merged.bam
#   04_pbulk/pseudobulk.cleaned.bam
#   04_pbulk/pseudobulk.bam
#   04_pbulk/pseudobulk.bam.bai
#   # 5.  Call variants on the pseudobulk sample.
#   05_call/interval_{interval}.intervals
#   05_call/variants.{interval}.vcf
#   05_call/variants.{interval}.snp_only.vcf
#   05_call/variants.table.txt
#   05_call/variants.metadata.txt
#   05_call/variants.unfiltered.matrix.{batch}.txt
#   05_call/variants.unfiltered.matrix.txt
#   05_call/variants.matrix.txt
#   05_call/variants.coord.txt
#   # 6.  Make a site x cell matrix containing the ref/alt read depth.
#   06_matrix/{sample}.{batch}.raw.pileup/
#   06_matrix/{sample}.{batch}.pileup/
#   06_matrix/{sample}.{batch}.coverage.vcf/
#   06_matrix/{sample}.{batch}.coverage.txt/
#   06_matrix/{sample}.{batch}.coverage.with_cells.txt
#   06_matrix/{sample}.{batch}.coverage.merged.txt
#   06_matrix/coverage.table.txt
#   06_matrix/coverage.matrix.txt
#   # 7.  Filter the sites.
#   07_filter/coverage.matrix.{batch}.txt
#   07_filter/coverage.matrix.filter1.{batch}.txt
#   07_filter/coverage.matrix.filter1.txt
#   07_filter/coverage.matrix.filter2.txt
#   07_filter/coverage.metadata.txt
#   07_filter/ref_count.all.txt
#   07_filter/alt_count.all.txt
#   07_filter/knn.features.txt
#   07_filter/ref_count.features.txt
#   07_filter/alt_count.features.txt
#   07_filter/knn.neighbors.{batch}.txt
#   07_filter/knn.neighbors.txt
#   07_filter/ref_count.{batch}.txt
#   07_filter/alt_count.{batch}.txt
#   07_filter/genotypes.{batch}.txt
#   07_filter/probabilities.{batch}.txt
#   07_filter/genotypes.txt
#   07_filter/probabilities.txt
#   07_filter/genotype.count.txt
#   07_filter/mixed_genotypes.txt
#   07_filter/coverage.matrix.filter3.txt
#   07_filter/coverage.matrix.filter4.txt
#   # 8.  Call/impute genotypes.
#   08_genotype/coverage.metadata.txt
#   08_genotype/ref_count.all.txt
#   08_genotype/ref_count.all.txt
#   08_genotype/knn.features.txt
#   08_genotype/ref_count.features.txt
#   08_genotype/alt_count.features.txt
#   08_genotype/knn.neighbors.0.txt
#   08_genotype/knn.neighbors.txt
#   08_genotype/genotypes.txt
#   08_genotype/probabilities.txt
#   # 9.  Make the phylogeny.
#   mutations.fa
#   phylogeny/
#     beast2.infile.txt
#     beast2.outfile.txt
#     beast2.model.RDS
#     screen.log
#     trace.log
#     tree.log
# logs/
#   genome.fa.fai.log
#   genome.dict.log
#   {sample}.alignments.log
#   {sample}.{batch}.alignments.log
#   {sample}.{batch}.cells.rg.log
#   {sample}.{batch}.cells.rg_contig.log
#   {sample}.{batch}.cells.rg_contig_index.log
#   {sample}.{batch}.cells.rg_contig_index_snt.log
#   {sample}.{batch}.cells.rg_contig_index_snt_report.log
#   {sample}.{batch}.cells.rg_contig_index_snt_recal.log
#   # 4.  Merge single cells to a pseudobulk sample.
#   {sample}.{batch}.pb.log
#   {sample}.pb.log
#   {sample}.pb.clean.log
#   pseudobulk.merged.log
#   pseudobulk.cleaned.log
#   pseudobulk.sort.log
#   pseudobulk.index.log
#   # 5.  Call variants on the pseudobulk sample.
#   interval_{interval}.intervals
#   pseudobulk.clean.contig.{interval}.log
# temp/
# scripts/
#   # 2.  Demultiplex CellRanger output into single cells.
#   extract_cellranger_bam.py
#   clean_cells_file.py
#   extract_sample_barcodes.py
#   extract_batch_barcodes.py
#   demux_one_batch.py
#   # 4.  Merge single cells to a pseudobulk sample.
#   clean_batch_header.py
#   # 5.  Call variants on the pseudobulk sample.
#   make_interval_file.py
#   keep_only_snps.py
#   parse_vcf_files.py
#   collect_variant_metadata.py
#   convert_to_matrix_by_batch.py
#   merge_matrix_batches.py
#   filter_variant_calls.py
#   extract_matrix_coordinates.py
#   # 6.  Make a site x cell matrix containing the ref/alt read depth.
#   extract_coverage_from_vcf.py
#   make_coverage_matrix.py
#   # 7.  Filter the sites.
#   split_matrix_by_row.py
#   filter_sites1.py
#   merge_filter1_matrix.py
#   filter_sites2.py
#   prepare_geno_calling_files1.py
#   select_geno_calling_features1.py
#   extract_geno_calling_features1.py
#   calc_neighbor_scores1.R
#   select_neighbors1.py
#   split_matrix_by_row.py
#   call_genotypes1.R
#   merge_matrix_by_row.py
#   count_genotypes.py
#   list_sites_with_mixed_genotypes.py
#   filter_sites3.py
#   filter_sites4.py
#   # 9.  Make the phylogeny.
#   make_fasta_file.py
#   run_beast2.R
#   summarize_beast2.R
#   analyze_mcc_tree.R
#   plot_model_estimates.R
#   plot_tree.R
#   plot_densitree.R



############################################################
# Rules
# -----
# You should not need to change anything below.

import os
opj = os.path.join


CELL_FILE = "data/cells.txt"
assert os.path.exists(CELL_FILE), "File not found: %s" % CELL_FILE

GENOME_DIR = "output"
DEMUX_DIR = "output/02_demux"
PREPROC_DIR = "output/03_preproc"
PBULK_DIR = "output/04_pbulk"
CALL_DIR = "output/05_call"
MATRIX_DIR = "output/06_matrix"
FILTER_DIR = "output/07_filter"
GENOTYPE_DIR = "output/08_genotype"
PHYLO_DIR = "output/phylogeny"

GENOME_LOG_DIR = "logs/01_genome"
DEMUX_LOG_DIR = "logs/02_demux"
PREPROC_LOG_DIR = "logs/03_preproc"
PBULK_LOG_DIR = "logs/04_pbulk"
CALL_LOG_DIR = "logs/05_call"
MATRIX_LOG_DIR = "logs/06_matrix"
PHYLO_LOG_DIR = "logs/phylogeny"


# Create the directories for the files we need to create.
x = ["output", GENOME_DIR, DEMUX_DIR, PREPROC_DIR, PBULK_DIR, CALL_DIR, 
     MATRIX_DIR, FILTER_DIR, GENOTYPE_DIR, PHYLO_DIR,
     "logs", GENOME_LOG_DIR, DEMUX_LOG_DIR, PREPROC_LOG_DIR, PBULK_LOG_DIR,
     CALL_LOG_DIR, MATRIX_LOG_DIR, PHYLO_LOG_DIR,
     "temp",
     ]
for x in x:
    if not os.path.exists(x):
        os.mkdir(x)




# Set default parameters, if not provided.

if BEAST2_DIR is None:
    # See if I can find the directory where BEAST2 is installed.
    # beast/bin/
    #   beast
    #   logcombiner
    #   treeannotator
    #   ...
    HOME = os.path.expanduser("~")
    BEAST2_DIR = os.path.join(
        HOME, "mambaforge/envs/snakemake/share/beast2-2.6.3-2")
    assert os.path.exists(BEAST2_DIR), "I could not find BEAST2_DIR"

if LOGCOMBINER is None:
    LOGCOMBINER = os.path.join(BEAST2_DIR, "bin", "logcombiner")
    assert os.path.exists(LOGCOMBINER), "I could not find LOGCOMBINER"
   



def read_cell_file(filename):
    """Read the user's cells.txt file to get a list of the samples.
    <sample>_<barcode>[-<gem_well>]   006_R_LUNG_CTGAATGCATCTTTCA
                                      006_R_LUNG_CTGAATGCATCTTTCA-1

    Return a dictionary of sample -> set of cells

    """
    i_cell = None
    cells = set()
    for line in open(filename):
        cols = line.rstrip("\r\n").split("\t")
        if i_cell is None:
            assert "Cell" in cols, \
                'File %s is missing a column with header "Cell".' % filename
            i_cell = cols.index("Cell")
            continue
        x = cols[i_cell].strip()
        if not x:
            continue
        cells.add(x)

    sample2cells = {}
    for cell in cells:
        x = cell.rsplit("_", 1)
        assert len(x) == 2, \
            "Cell not formatted as <sample>_<barcode>: %s" % cell
        sample2cells.setdefault(x[0], set()).add(cell)

    return sample2cells



# CELLS_PER_BATCH_PREPROCESSING

PREPROCESSING_BATCHES = []   # list of (sample, batch)
sample2cells = read_cell_file(CELL_FILE)
for sample, cells in sample2cells.items():
    num_batches = int((len(cells)-1)/CELLS_PER_BATCH_PREPROCESSING) + 1
    nd = len(str(num_batches))
    for i in range(num_batches):
        batch = "%0*d" % (nd, i)
        PREPROCESSING_BATCHES.append((sample, batch))
x = [len(x) for x in sample2cells.values()]
NUM_CELLS = sum(x)







# CELLS_PER_GENOME_INTERVAL

x = NUM_CELLS / CELLS_PER_GENOME_INTERVAL
n = max(int(x), 1)
nd = len(str(n))
INTERVALS = ["%0*d" % (nd, i) for i in range(n)]



# CELLS_PER_BATCH_VARIANT_MATRIX
# We want to do one batch per 64 Mb of pseudobulk.variants.txt
# file.
# 14 cells generated a 34 Mb file.

x = NUM_CELLS / CELLS_PER_BATCH_VARIANT_MATRIX
n = max(int(x), 1)
nd = len(str(n))
MATRIX_BATCHES = ["%0*d" % (nd, i) for i in range(n)]



# CELLS_PER_BATCH_FILTER1
# Want to filter 256 Mb at a time.
# 14 cells generated a 5.1 Mb file.

x = NUM_CELLS / CELLS_PER_BATCH_FILTER1
n = max(int(x), 1)
nd = len(str(n))
FILTER1_BATCHES = ["%0*d" % (nd, i) for i in range(n)]



# CELLS_PER_BATCH_KNN
# Use ~512 cells / batch.
# Not sure exactly how many cells here, due to filtering.

x = NUM_CELLS / CELLS_PER_BATCH_KNN
n = max(int(x), 1)
nd = len(str(n))
KNN_BATCHES = ["%0*d" % (nd, i) for i in range(n)]



# CELLS_PER_BATCH_GENOTYPE_CALLING
# Want ~64 Mb / batch.
# 14 cells generated a 2.3 Kb file, so extrapolate, around 390 cells.

x = NUM_CELLS / CELLS_PER_BATCH_GENOTYPE_CALLING
n = max(int(x), 1)
nd = len(str(n))
GENOTYPE_CALLING_BATCHES = ["%0*d" % (nd, i) for i in range(n)]




rule all:
    input:
        opj(PHYLO_DIR, "max_clade_cred.nexus.txt"),
        opj(PHYLO_DIR, "max_clade_cred.newick.txt"),
        opj(PHYLO_DIR, "max_clade_cred.metadata.txt"),
        opj(PHYLO_DIR, "summary.txt"),
        opj(PHYLO_DIR, "posterior.pdf"),
        opj(PHYLO_DIR, "tree.mcc.pdf"),
        opj(PHYLO_DIR, "densitree.pdf"),


rule copy_ref_genome:
    input:
        "data/genome.fa"
    output:
        opj(GENOME_DIR, "genome.fa")
    params:
        GENOME_DIR=GENOME_DIR,
    shell:
        "mkdir -p {params.GENOME_DIR}; cp -p {input} {output}"


rule index_ref_genome:
    input:
        opj(GENOME_DIR, "genome.fa")
    output:
        opj(GENOME_DIR, "genome.fa.fai")
    log:
        opj(GENOME_LOG_DIR, "genome.fa.fai.log")
    params:
        SAMTOOLS=SAMTOOLS,
    shell:
        "{params.SAMTOOLS} faidx {input} --fai-idx {output} >& {log}"


rule create_ref_genome_dict:
    input:
        opj(GENOME_DIR, "genome.fa")
    output:
        opj(GENOME_DIR, "genome.dict")
    log:
        opj(GENOME_LOG_DIR, "genome.dict.log")
    params:
        PICARD=PICARD,
    shell:
        """
        {params.PICARD} CreateSequenceDictionary \
             R={input} \
             O={output} >& {log}
        """


rule extract_cellranger_bam:
    input:
        "data/cellranger"
    output:
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.bam")
    params:
        sample="{sample}"
    script:
        "scripts/extract_cellranger_bam.py"


rule extract_bam_header:
    input:
        opj(DEMUX_DIR, "{sample}.bam")
    output:
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.header.txt")
    params:
        SAMTOOLS=SAMTOOLS
    shell:
        "{params.SAMTOOLS} view {input} -H > {output}"


rule clean_cells_file:
    input:
        CELL_FILE,
        list(opj(DEMUX_DIR, "%s.bam" % sample)
             for (sample, batch) in PREPROCESSING_BATCHES),
    output:
        opj(DEMUX_DIR, "cells.txt")
    script:
        "scripts/clean_cells_file.py"


rule extract_sample_barcodes:
    input:
        opj(DEMUX_DIR, "cells.txt")
    output:
        # Do not allow periods in the name to disambiguate from
        # extract_barcodes_for_batch.
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.barcodes.txt"),
    params:
        sample="{sample}",
    script:
        "scripts/extract_sample_barcodes.py"


rule extract_sample_alignments:
    input:
        bam_file=opj(DEMUX_DIR, "{sample}.bam"),
        barcode_file=opj(DEMUX_DIR, "{sample}.barcodes.txt"),
    output:
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.alignments.txt")
    log:
        opj(DEMUX_LOG_DIR, "{sample}.alignments.log")
    params:
        SAMTOOLS=SAMTOOLS
    shell:
         "{params.SAMTOOLS} view {input.bam_file} | \
          LC_ALL=C grep -F -f {input.barcode_file} 1> {output} 2> {log}"


rule extract_batch_barcodes:
    input:
        opj(DEMUX_DIR, "cells.txt")
    output:
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.barcodes.txt"),
    params:
        sample="{sample}",
        batch="{batch}",
        batch_size=CELLS_PER_BATCH_PREPROCESSING,
    script:
        "scripts/extract_batch_barcodes.py"


rule extract_batch_alignments:
    input:
        align_file=opj(DEMUX_DIR, "{sample}.alignments.txt"),
        barcode_file=opj(DEMUX_DIR, "{sample}.{batch}.barcodes.txt"),
    output:
        opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.alignments.txt")
    log:
        opj(DEMUX_LOG_DIR, "{sample}.{batch}.alignments.log")
    shell:
        "cat {input.align_file} | \
         LC_ALL=C grep -F -f {input.barcode_file} 1> {output} 2> {log}"


rule demux_one_batch:
    input:
        align_file=opj(DEMUX_DIR, "{sample}.{batch}.alignments.txt"),
        header_file=opj(DEMUX_DIR, "{sample}.header.txt"),
    output:
        directory(
            opj(DEMUX_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.sam")),
    params:
        sample="{sample}",
    script:
        "scripts/demux_one_batch.py"


rule convert_sam_to_bam:
    input:
        opj(DEMUX_DIR, "{sample}.{batch}.cells.sam")
    output:
        directory(
            opj(PREPROC_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.bam"))
    params:
        DEMUX_DIR=DEMUX_DIR,
        PREPROC_DIR=PREPROC_DIR,
        SAMTOOLS=SAMTOOLS,
    shell:
        """mkdir -p {params.PREPROC_DIR}
        for i in {input}/*.sam; do
             j=`echo $i | sed -e 's/.sam/.bam/g'`
             j=`echo $j | sed -e 's#{params.DEMUX_DIR}#{params.PREPROC_DIR}#g'`
             mkdir -p `dirname $j`
             {params.SAMTOOLS} view -bS -o $j $i
        done
        """


rule add_read_groups:
    input:
        opj(PREPROC_DIR, "{sample}.{batch}.cells.bam")
    output:
        directory(opj(
            PREPROC_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg.bam"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg.log")
    params:
        #sample="{sample}",
        PICARD=PICARD,
    shell:
        """
        for i in {input}/*.bam; do 
             j=`echo $i | sed -e 's/cells.bam/cells.rg.bam/'`
             mkdir -p `dirname $j`
             {params.PICARD} AddOrReplaceReadGroups \
                 I=$i \
                 O=$j \
                 ID=group \
                 LB=library \
                 PU=platform \
                 SM=pseudobulk \
                 PL=ILLUMINA \
                 VALIDATION_STRINGENCY=LENIENT >& {log}
        done
        """


rule sort_bam:
    input:
        bam_dir=opj(PREPROC_DIR, "{sample}.{batch}.cells.rg.bam"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
    output:
        directory(opj(
            PREPROC_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg_contig.bam"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg_contig.log")
    params:
        PICARD=PICARD
    shell:
        # Picard version >= 2.24 requires REFERENCE_SEQUENCE= and
        # SEQUENCE_DICTIONARY=.
        # On 230606, version in BioConda is 2.18.29-SNAPSHOT.
        # On 230721, version in BioConda is 3.0.0.
        """
        for i in {input.bam_dir}/*.bam; do
             j=`echo $i | sed -e 's/cells.rg.bam/cells.rg_contig.bam/'`
             mkdir -p `dirname $j`
             {params.PICARD} ReorderSam \
                 -I $i \
                 -O $j \
                 -R {input.ref_file} \
                 -SD {input.ref_dict} \
                 --VALIDATION_STRINGENCY LENIENT \
                 --ALLOW_INCOMPLETE_DICT_CONCORDANCE true \
                 --TMP_DIR temp >& {log}
        done
        # For some reason, Picard makes the permission of the TEMP_DIR
        # to 777.  Fix it.
        chmod 0755 temp
        """


rule index_bam:
    input:
        opj(PREPROC_DIR, "{sample}.{batch}.cells.rg_contig.bam")
    output:
        directory(opj(PREPROC_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg_contig_index.bam"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg_contig_index.log")
    params:
        SAMTOOLS=SAMTOOLS
    shell:
        """
        for i in {input}/*.bam; do
             j=`echo $i | sed -e 's/cells.rg_contig.bam/cells.rg_contig_index.bam/'`
             mkdir -p `dirname $j`
             cp -p $i $j
             {params.SAMTOOLS} index $j >& {log}
        done
        """


rule split_n_trim:
    input:
        bam_dir=opj(PREPROC_DIR, "{sample}.{batch}.cells.rg_contig_index.bam"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
    output:
        directory(opj(PREPROC_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg_contig_index_snt.bam"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg_contig_index_snt.log")
    params:
        GATK=GATK
    shell:
        """
        for i in {input.bam_dir}/*.bam; do
             j=`echo $i | sed -e 's/cells.rg_contig_index.bam/cells.rg_contig_index_snt.bam/'`
             mkdir -p `dirname $j`
             {params.GATK} SplitNCigarReads \
                 --max-reads-in-memory 150000 \
                 --tmp-dir temp \
                 -R {input.ref_file} \
                 -I $i \
                 -O $j >& {log}
        done"""


rule make_base_recalibration_report:
    input:
        bam_dir=opj(
            PREPROC_DIR, "{sample}.{batch}.cells.rg_contig_index_snt.bam"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
        # Make these an input rather than a parameter so
        # Snakemake will complain if the file is missing.
        recal_known_sites1="data/known_sites1.vcf.gz",
        recal_known_sites1_tbi="data/known_sites1.vcf.gz.tbi",
        recal_known_sites2="data/known_sites2.vcf.gz",
        recal_known_sites2_tbi="data/known_sites2.vcf.gz.tbi",
        recal_known_sites3="data/known_sites3.vcf.gz",
        recal_known_sites3_tbi="data/known_sites3.vcf.gz.tbi",
    output:
        directory(opj(PREPROC_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg_contig_index_snt.report"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg_contig_index_snt_report.log")
    params:
        GATK=GATK
    shell:
        """
        for i in {input.bam_dir}/*.bam; do
             j=`echo $i | sed -e 's/cells.rg_contig_index_snt.bam/cells.rg_contig_index_snt.report/'`
             mkdir -p `dirname $j`
             {params.GATK} BaseRecalibrator \
                 --reference {input.ref_file} \
                 --input $i \
                 --output $j \
                 --known-sites {input.recal_known_sites1} \
                 --known-sites {input.recal_known_sites2} \
                 --known-sites {input.recal_known_sites3} \
                 >& {log}
        done"""


rule recalibrate_base_quality_score:
    input:
        bam_dir=opj(
            PREPROC_DIR, "{sample}.{batch}.cells.rg_contig_index_snt.bam"),
        report_dir=opj(
            PREPROC_DIR, "{sample}.{batch}.cells.rg_contig_index_snt.report"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
    output:
        directory(opj(PREPROC_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.cells.rg_contig_index_snt_recal.bam"))
    log:
        opj(PREPROC_LOG_DIR, "{sample}.{batch}.cells.rg_contig_index_snt_recal.log")
    params:
        GATK=GATK
    shell:
        """
        for i in {input.bam_dir}/*.bam; do
            j=`echo $i | sed -e 's/cells.rg_contig_index_snt.bam/cells.rg_contig_index_snt_recal.bam/'`
            k=`echo $i | sed -e 's/cells.rg_contig_index_snt.bam/cells.rg_contig_index_snt.report/'`
            mkdir -p `dirname $j`
            {params.GATK} ApplyBQSR \
                -R {input.ref_file} \
                --bqsr-recal-file $k \
                -I $i \
                -O $j \
                >& {log}
        done"""


rule merge_cells_to_batch:
    input:
        opj(PREPROC_DIR,
            "{sample}.{batch}.cells.rg_contig_index_snt_recal.bam")
    output:
        opj(PBULK_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.merged.bam")
    log:
        opj(PBULK_LOG_DIR, "{sample}.{batch}.merged.log")
    params:
        sample="{sample}",
        batch="{batch}",
        PBULK_DIR=PBULK_DIR,
        SAMTOOLS=SAMTOOLS,
    shell:
        """
        mkdir -p {params.PBULK_DIR}
        FILES=`ls {input}/*.bam`
        TF={params.PBULK_DIR}/{params.sample}.{params.batch}.files
        echo "${{FILES[*]}}" > $TF
        {params.SAMTOOLS} merge -f -b $TF {output} >& {log}
        """


rule clean_batch_header:
    input:
        opj(PBULK_DIR, "{sample}.{batch}.merged.bam")
    output:
        opj(PBULK_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.bam")
    params:
        header_file=opj(PBULK_DIR, "{sample}.{batch}.merged.header"),
        clean_header_file=opj(PBULK_DIR, "{sample}.{batch}.header"),
        log_file=opj(PBULK_LOG_DIR, "{sample}.{batch}.pb.clean.log"),
    script:
        "scripts/clean_batch_header.py"


rule merge_batches_to_sample:
    input:
        lambda wildcards:
            [opj(PBULK_DIR, "{sample}.%s.bam" % batch)
             for (sample, batch) in PREPROCESSING_BATCHES
             if sample == wildcards.sample]
    output:
        opj(PBULK_DIR, "{sample,[A-Za-z0-9_-]+}.merged.bam")
    log:
        opj(PBULK_LOG_DIR, "{sample}.merged.log")
    params:
        sample="{sample}",
        SAMTOOLS=SAMTOOLS,
        PBULK_DIR=PBULK_DIR,
    shell: 
        """
        TF={params.PBULK_DIR}/{params.sample}.files
        echo "{input}" | tr ' ' '\n' > $TF
        {params.SAMTOOLS} merge -f -b $TF {output} >& {log}
        """


rule clean_sample_header:
    input:
        opj(PBULK_DIR, "{sample}.merged.bam")
    output:
        opj(PBULK_DIR, "{sample,[A-Za-z0-9_-]+}.bam")
    params:
        header_file=opj(PBULK_DIR, "{sample}.merged.header"),
        clean_header_file=opj(PBULK_DIR, "{sample}.header"),
        log_file=opj(PBULK_LOG_DIR, "{sample}.log"),
    script:
        "scripts/clean_batch_header.py"


rule merge_samples_to_pseudobulk:
    input:
        list(opj(PBULK_DIR, "%s.bam" % sample)
            for (sample, batch) in PREPROCESSING_BATCHES)
    output:
        opj(PBULK_DIR, "pseudobulk.merged.bam"),
    log:
        opj(PBULK_LOG_DIR, "pseudobulk.merged.log")
    params:
        SAMTOOLS=SAMTOOLS,
        PBULK_DIR=PBULK_DIR,
    shell: 
        """
        TF={params.PBULK_DIR}/pseudobulk.files
        echo "{input}" | tr ' ' '\n' > $TF
        {params.SAMTOOLS} merge -f -b $TF {output} >& {log}
        """

rule clean_pseudobulk_header:
    input:
        opj(PBULK_DIR, "pseudobulk.merged.bam"),
    output:
        opj(PBULK_DIR, "pseudobulk.cleaned.bam"),
    params:
        header_file=opj(PBULK_DIR, "pseudobulk.merged.header"),
        clean_header_file=opj(PBULK_DIR, "pseudobulk.header"),
        log_file=opj(PBULK_LOG_DIR, "pseudobulk.cleaned.log"),
    script:
        "scripts/clean_batch_header.py"


rule sort_pseudobulk_by_contig:
    input:
        bam_file=opj(PBULK_DIR, "pseudobulk.cleaned.bam"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
    output:
        opj(PBULK_DIR, "pseudobulk.bam")
    log:
        opj(PBULK_LOG_DIR, "pseudobulk.sort.log")
    params:
        PICARD=PICARD
    shell:
        # Picard version >= 2.24 requires REFERENCE_SEQUENCE= and
        # SEQUENCE_DICTIONARY=.
        # Current version (on 230606) in BioConda is
        # 2.18.29-SNAPSHOT.
        """
        {params.PICARD} -Xmx16g ReorderSam \
             -I {input.bam_file} \
             -O {output} \
             -R {input.ref_file} \
             -SD {input.ref_dict} \
             --VALIDATION_STRINGENCY LENIENT \
             --ALLOW_INCOMPLETE_DICT_CONCORDANCE true \
             --TMP_DIR temp >& {log}
        # For some reason, Picard makes the permission of the TEMP_DIR
        # to 777.  Fix it.
        chmod 0755 temp
        """


rule index_pseudobulk_bam:
    input:
        opj(PBULK_DIR, "pseudobulk.bam")
    output:
        opj(PBULK_DIR, "pseudobulk.bam.bai")
    log:
        opj(PBULK_LOG_DIR, "pseudobulk.index.log")
    params:
        SAMTOOLS=SAMTOOLS
    shell:
        "{params.SAMTOOLS} index {input} >& {log}"


rule make_interval_file:
    input:
        opj(GENOME_DIR, "genome.fa.fai")
    output:
        opj(CALL_DIR, "interval_{interval,\d+}.intervals")
    params:
        interval="{interval}",
        num_intervals=len(INTERVALS),
    script:
        "scripts/make_interval_file.py"


rule call_variants:
    input:
        bam_file=opj(PBULK_DIR, "pseudobulk.bam"),
        bam_index=opj(PBULK_DIR, "pseudobulk.bam.bai"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
        interval_file=opj(CALL_DIR, "interval_{interval}.intervals"),
    output:
        opj(CALL_DIR, "variants.{interval,\d+}.vcf")
    log:
        opj(CALL_LOG_DIR, "variants.{interval}.log")
    params:
        GATK=GATK
    shell:
        """
        {params.GATK} HaplotypeCaller \
            --input {input.bam_file} \
            --output {output} \
            --reference {input.ref_file} \
            --intervals {input.interval_file} \
            --interval-padding 150 \
            --dont-use-soft-clipped-bases true >& {log}
        """


rule keep_only_snps:
    input:
        opj(CALL_DIR, "variants.{interval}.vcf")
    output:
        opj(CALL_DIR, "variants.{interval,\d+}.snp_only.vcf")
    script:
        "scripts/keep_only_snps.py"


rule parse_vcf_files:
    input:
        [opj(CALL_DIR, "variants.%s.snp_only.vcf" % x) for x in INTERVALS]
    output:
        opj(CALL_DIR, "variants.table.txt")
    script:
        "scripts/parse_vcf_files.py"


rule collect_variant_metadata:
    input:
        opj(CALL_DIR, "variants.table.txt")
    output:
        opj(CALL_DIR, "variants.metadata.txt")
    script:
        "scripts/collect_variant_metadata.py"


rule convert_to_matrix_by_batch:
    input:
        opj(CALL_DIR, "variants.table.txt"),
        opj(CALL_DIR, "variants.metadata.txt"),
    output:
        opj(CALL_DIR, "variants.unfiltered.matrix.{batch,[A-Za-z0-9_-]+}.txt")
    params:
        batch="{batch}",
        num_batches=len(MATRIX_BATCHES),
    script:
        "scripts/convert_to_matrix_by_batch.py"


rule merge_matrix_batches:
    input:
        [opj(CALL_DIR, "variants.unfiltered.matrix.%s.txt" % x)
         for x in MATRIX_BATCHES]
    output:
        opj(CALL_DIR, "variants.unfiltered.matrix.txt")
    script:
        "scripts/merge_matrix_batches.py"


rule filter_variant_calls:
    input:
        opj(CALL_DIR, "variants.unfiltered.matrix.txt")
    output:
        opj(CALL_DIR, "variants.matrix.txt")
    params:
        keep_calls_with_min_alt_reads=KEEP_CALLS_WITH_MIN_ALT_READS,
        keep_calls_with_min_total_reads=KEEP_CALLS_WITH_MIN_TOTAL_READS,
        keep_calls_with_min_vaf=KEEP_CALLS_WITH_MIN_VAF,
    script:
        "scripts/filter_variant_calls.py"


rule extract_matrix_coordinates:
    input:
        opj(CALL_DIR, "variants.matrix.txt")
    output:
        opj(CALL_DIR, "variants.coord.txt")
    script:
        "scripts/extract_matrix_coordinates.py"


rule make_pileup_from_cells:
    input:
        bam_dir=opj(PREPROC_DIR,
            "{sample}.{batch}.cells.rg_contig_index_snt_recal.bam"),
        coord_file=opj(CALL_DIR, "variants.coord.txt"),
        ref_file=opj(GENOME_DIR, "genome.fa"),
        ref_index=opj(GENOME_DIR, "genome.fa.fai"),
        ref_dict=opj(GENOME_DIR, "genome.dict"),
    output:
        directory(opj(
            MATRIX_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.raw.pileup"))
    log:
        opj(MATRIX_LOG_DIR, "{sample}.{batch}.raw.pileup.log")
    params:
        PREPROC_DIR=PREPROC_DIR,
        MATRIX_DIR=MATRIX_DIR,
        SAMTOOLS=SAMTOOLS,
    shell:
        """
        for i in {input.bam_dir}/*.bam; do
             j=`echo $i | \
               sed -e 's/.cells.rg_contig_index_snt_recal.bam/.raw.pileup/'`
             j=`echo $j | sed -e 's/\.bam$/.pileup/'`
             j=`echo $j | \
               sed -e 's#{params.PREPROC_DIR}#{params.MATRIX_DIR}#g'`
             k=`echo $j | sed -e 's#output/#logs/#'`
             mkdir -p `dirname $j`
             mkdir -p `dirname $k`
             {params.SAMTOOLS} mpileup \
                  -f {input.ref_file} \
                  -l {input.coord_file} \
                  -R -B -q 0 -Q 0 -d10000000 \
                  $i 2> $k 1> $j
        done >& {log}
        """



rule clean_pileup_file:
    # VarScan will generate a "Parsing Exception" if there are 0
    # reads in a location.  Filter those lines out.
    input:
        opj(MATRIX_DIR, "{sample}.{batch}.raw.pileup"),
    output:
        directory(opj(
            MATRIX_DIR, "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.pileup")),
    shell:
        """
        for i in {input}/*.pileup; do 
             j=`echo $i | sed -e 's/.raw.pileup/.pileup/'`
             mkdir -p `dirname $j`
             cat $i | awk -F'\t' '$4 != 0 {{print}}' > $j
        done
        """


rule convert_pileup_to_vcf:
    input:
        opj(MATRIX_DIR, "{sample}.{batch}.pileup")
    output:
        directory(opj(MATRIX_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.coverage.vcf"))
    params:
        VARSCAN=VARSCAN
    shell:
        """
        for i in {input}/*.pileup; do
             j=`echo $i | sed -e 's/.pileup/.coverage.vcf/'`
             j=`echo $j | sed -e 's/.pileup$/.vcf/'`
             k=`echo $j | sed -e 's#output/#logs/#'`
             k=`echo $k | sed -e 's/.vcf$/.log/'`
             mkdir -p `dirname $j`
             mkdir -p `dirname $k`
             {params.VARSCAN} mpileup2cns $i \
                --min-coverage 0 \
                --min-reads2 0 \
                --min-avg-qual 0 \
                --min-var-freq 0 \
                --p-value 1.0 \
                --strand-filter 0 \
                --output-vcf 1 \
                1> $j \
                2> $k
        done
        """


rule extract_coverage_from_vcf:
    input:
        opj(MATRIX_DIR, "{sample}.{batch}.coverage.vcf")
    output:
        directory(opj(MATRIX_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.coverage.txt"))
    params:
        sample="{sample}",
    script:
        "scripts/extract_coverage_from_vcf.py"


rule add_cell_names_to_coverage:
    input:
        opj(MATRIX_DIR, "{sample}.{batch}.coverage.txt")
    output:
        directory(opj(MATRIX_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.coverage.with_cells.txt"))
    shell:
        # Input looks like:
        # output/015_R_LIVER.0.coverage.txt/015_R_LIVER_TGACAGTAGTGACCTT-1.txt
        # Pull out the cell name:
        # 015_R_LIVER_TGACAGTAGTGACCTT-1
        """
        for i in {input}/*.txt; do
            j=`echo $i | sed -e 's/.coverage.txt/.coverage.with_cells.txt/'`
            k=`basename $i`
            k=`echo $k | sed -e 's/.txt$//'`
            mkdir -p `dirname $j`
            sed "1s/\$/\tCell/; 2,\$s/\$/\t$k/" $i > $j
        done
        """


rule merge_coverage_by_batch:
    input:
        opj(MATRIX_DIR, "{sample}.{batch}.coverage.with_cells.txt")
    output:
        opj(MATRIX_DIR,
            "{sample,[A-Za-z0-9_-]+}.{batch,\d+}.coverage.merged.txt")
    log:
        opj(MATRIX_LOG_DIR, "{sample}.{batch}.coverage.merged.log")
    params:
        MATRIX_DIR=MATRIX_DIR,
        sample="{sample}",
        batch="{batch}",
    shell:
        """
        temp1={params.MATRIX_DIR}/{params.sample}.{params.batch}.coverage.merged.01.txt
        temp2={params.MATRIX_DIR}/{params.sample}.{params.batch}.coverage.merged.02.txt
        temp3={params.MATRIX_DIR}/{params.sample}.{params.batch}.coverage.merged.03.txt
        cat {input}/*.txt > $temp1
        grep -v "Chrom" $temp1 > $temp2
        sort -T . --parallel 4 -k1,1 -k2,2n -k3,4 $temp2 1> $temp3 2> {log}
        cat <(head -1 $temp1) $temp3 > {output}
        """


rule merge_coverage_batches:
    input:
        [opj(MATRIX_DIR, "%s.%s.coverage.merged.txt" % (sample, batch))
         for (sample, batch) in PREPROCESSING_BATCHES]
    output:
        opj(MATRIX_DIR, "coverage.table.txt")
    params:
        MATRIX_DIR=MATRIX_DIR
    shell:
        """
        temp1={params.MATRIX_DIR}/coverage.01.txt
        temp2={params.MATRIX_DIR}/coverage.02.txt
        temp3={params.MATRIX_DIR}/coverage.03.txt
        cat {input} > $temp1
        grep -v "Chrom" $temp1 > $temp2
        sort -T . --parallel 4 -k1,1 -k2,2n -k3,4 $temp2 > $temp3
        cat <(head -1 $temp1) $temp3 > {output}
        """


rule make_coverage_matrix:
    input:
        opj(MATRIX_DIR, "coverage.table.txt"),
        opj(CALL_DIR, "variants.matrix.txt"),
    output:
        opj(MATRIX_DIR, "coverage.matrix.txt")
    script:
        "scripts/make_coverage_matrix.py"


rule split_matrix_for_filter1:
    input:
        opj(MATRIX_DIR, "coverage.matrix.txt"),
    output:
        [opj(FILTER_DIR, "coverage.matrix.%s.txt" % batch)
         for batch in FILTER1_BATCHES]
    script:
        "scripts/split_matrix_by_row.py"


rule filter_sites1:
    input:
        opj(FILTER_DIR, "coverage.matrix.{batch,[A-Za-z0-9_-]+}.txt")
    output:
        opj(FILTER_DIR, "coverage.matrix.filter1.{batch,[A-Za-z0-9_-]+}.txt")
    params:
        keep_sites_with_high_reads=KEEP_SITES_WITH_HIGH_READS,
        remove_chrom_with_prefix=REMOVE_CHROM_WITH_PREFIX,
        remove_multiallelic_sites=REMOVE_MULTIALLELIC_SITES,
    script:
        "scripts/filter_sites1.py"


rule merge_filter1_matrix:
    input:
        [opj(FILTER_DIR, "coverage.matrix.filter1.%s.txt" % batch)
         for batch in FILTER1_BATCHES],
    output:
        opj(FILTER_DIR, "coverage.matrix.filter1.txt")
    script:
        "scripts/merge_filter1_matrix.py"


rule filter_sites2:
    input:
        opj(FILTER_DIR, "coverage.matrix.filter1.txt")
    output:
        opj(FILTER_DIR, "coverage.matrix.filter2.txt")
    params:
        remove_close_sites=REMOVE_CLOSE_SITES,
    script:
        "scripts/filter_sites2.py"


rule prepare_geno_calling_files1:
    input:
        opj(FILTER_DIR, "coverage.matrix.filter2.txt")
    output:
        opj(FILTER_DIR, "coverage.metadata.txt"),
        opj(FILTER_DIR, "ref_count.all.txt"),
        opj(FILTER_DIR, "alt_count.all.txt"),
    script:
        "scripts/prepare_geno_calling_files1.py"


rule select_geno_calling_features1:
    input:
        opj(FILTER_DIR, "coverage.metadata.txt"),
    output:
        opj(FILTER_DIR, "knn.features.txt"),
    params:
        ternary_features=TERNARY_FEATURES,
    script:
        "scripts/select_geno_calling_features1.py"


rule extract_geno_calling_features1:
    input:
        opj(FILTER_DIR, "ref_count.all.txt"),
        opj(FILTER_DIR, "alt_count.all.txt"),
        opj(FILTER_DIR, "knn.features.txt"),
    output:
        opj(FILTER_DIR, "ref_count.features.txt"),
        opj(FILTER_DIR, "alt_count.features.txt"),
    script:
        "scripts/extract_geno_calling_features1.py"


rule calc_neighbor_scores1:
    input:
        opj(FILTER_DIR, "ref_count.features.txt"),
        opj(FILTER_DIR, "alt_count.features.txt"),
    output:
        opj(FILTER_DIR, "knn.neighbors.{batch,[A-Za-z0-9_-]+}.txt"),
    params:
        K=TERNARY_K,
        num_samples=TERNARY_SAMPLES,
        rng_seed=PHYLO_RNG_SEED,
        start=int(list({batch})[0])+1,
        skip=len(KNN_BATCHES),
        num_cores=workflow.cores,
    script:
        "scripts/calc_neighbor_scores1.R"


rule select_neighbors1:
    input:
        [opj(FILTER_DIR, "knn.neighbors.%s.txt" % batch)
         for batch in KNN_BATCHES],
    output:
        opj(FILTER_DIR, "knn.neighbors.txt"),
    params:
        K=TERNARY_K,
    script:
        "scripts/select_neighbors1.py"


rule split_ref_matrix1:
    input:
        opj(FILTER_DIR, "ref_count.all.txt"),
    output:
        [opj(FILTER_DIR, "ref_count.%s.txt" % batch)
         for batch in GENOTYPE_CALLING_BATCHES],
    script:
        "scripts/split_matrix_by_row.py"


rule split_alt_matrix1:
    input:
        opj(FILTER_DIR, "alt_count.all.txt"),
    output:
        [opj(FILTER_DIR, "alt_count.%s.txt" % batch)
         for batch in GENOTYPE_CALLING_BATCHES],
    script:
        "scripts/split_matrix_by_row.py"


rule call_genotypes1:
    input:
        opj(FILTER_DIR, "ref_count.{batch}.txt"),
        opj(FILTER_DIR, "alt_count.{batch}.txt"),
        opj(FILTER_DIR, "knn.neighbors.txt"),
    output:
        opj(FILTER_DIR, "genotypes.{batch,[A-Za-z0-9_-]+}.txt"),
        opj(FILTER_DIR, "probabilities.{batch,[A-Za-z0-9_-]+}.txt"),
    params:
        delta=TERNARY_DELTA,
        impute=TERNARY_IMPUTE,
        num_cores=workflow.cores,
    script:
        "scripts/call_genotypes1.R"


rule merge_genotype_files1:
    input:
        [opj(FILTER_DIR, "genotypes.%s.txt" % batch)
         for batch in GENOTYPE_CALLING_BATCHES],
    output:
        opj(FILTER_DIR, "genotypes.txt"),
    script:
        "scripts/merge_matrix_by_row.py"


rule merge_probabilities_files1:
    input:
        [opj(FILTER_DIR, "probabilities.%s.txt" % batch)
         for batch in GENOTYPE_CALLING_BATCHES]
    output:
        opj(FILTER_DIR, "probabilities.txt")
    script:
        "scripts/merge_matrix_by_row.py"


rule count_genotypes:
    input:
        opj(FILTER_DIR, "genotypes.txt")
    output:
        opj(FILTER_DIR, "genotype.count.txt"),
    script:
        "scripts/count_genotypes.py"


rule list_sites_with_mixed_genotypes:
    input:
        opj(FILTER_DIR, "genotype.count.txt"),
    output:
        opj(FILTER_DIR, "mixed_genotypes.txt"),
    params:
        keep_sites_with_mixed_calls=KEEP_SITES_WITH_MIXED_GENOTYPES,
    script:
        "scripts/list_sites_with_mixed_genotypes.py"


rule filter_sites3:
    input:
        opj(FILTER_DIR, "coverage.matrix.filter2.txt"),
        opj(FILTER_DIR, "mixed_genotypes.txt"),
    output:
        opj(FILTER_DIR, "coverage.matrix.filter3.txt")
    script:
        "scripts/filter_sites3.py"


rule filter_sites4:
    input:
        opj(FILTER_DIR, "coverage.matrix.filter3.txt")
    output:
        opj(FILTER_DIR, "coverage.matrix.filter4.txt")
    params:
        keep_n_sites_seen_in_most_cells=KEEP_N_SITES_SEEN_IN_MOST_CELLS,
    script:
        "scripts/filter_sites4.py"


rule prepare_geno_calling_files2:
    input:
        opj(FILTER_DIR, "coverage.matrix.filter4.txt")
    output:
        opj(GENOTYPE_DIR, "coverage.metadata.txt"),
        opj(GENOTYPE_DIR, "ref_count.all.txt"),
        opj(GENOTYPE_DIR, "alt_count.all.txt"),
    script:
        "scripts/prepare_geno_calling_files1.py"


rule select_geno_calling_features2:
    input:
        opj(GENOTYPE_DIR, "coverage.metadata.txt"),
    output:
        opj(GENOTYPE_DIR, "knn.features.txt"),
    params:
        ternary_features=TERNARY_FEATURES,
    script:
        "scripts/select_geno_calling_features1.py"


rule extract_geno_calling_features2:
    input:
        opj(GENOTYPE_DIR, "ref_count.all.txt"),
        opj(GENOTYPE_DIR, "alt_count.all.txt"),
        opj(GENOTYPE_DIR, "knn.features.txt"),
    output:
        opj(GENOTYPE_DIR, "ref_count.features.txt"),
        opj(GENOTYPE_DIR, "alt_count.features.txt"),
    script:
        "scripts/extract_geno_calling_features1.py"


rule calc_neighbor_scores2:
    input:
        opj(GENOTYPE_DIR, "ref_count.features.txt"),
        opj(GENOTYPE_DIR, "alt_count.features.txt"),
    output:
        opj(GENOTYPE_DIR, "knn.neighbors.0.txt"),
    params:
        K=TERNARY_K,
        num_samples=TERNARY_SAMPLES,
        rng_seed=PHYLO_RNG_SEED,
        start=1,
        skip=1,
        num_cores=workflow.cores,
    script:
        "scripts/calc_neighbor_scores1.R"


rule select_neighbors2:
    input:
        opj(GENOTYPE_DIR, "knn.neighbors.0.txt"),
    output:
        opj(GENOTYPE_DIR, "knn.neighbors.txt"),
    params:
        K=TERNARY_K,
    script:
        "scripts/select_neighbors1.py"


rule call_genotypes2:
    input:
        opj(GENOTYPE_DIR, "ref_count.all.txt"),
        opj(GENOTYPE_DIR, "alt_count.all.txt"),
        opj(GENOTYPE_DIR, "knn.neighbors.txt"),
    output:
        opj(GENOTYPE_DIR, "genotypes.txt"),
        opj(GENOTYPE_DIR, "probabilities.txt"),
    params:
        delta=TERNARY_DELTA,
        impute=TERNARY_IMPUTE,
        num_cores=workflow.cores,
    script:
        "scripts/call_genotypes1.R"


rule make_fasta_file:
    input:
        opj(GENOTYPE_DIR, "genotypes.txt"),
    output:
        opj(PHYLO_DIR, "mutations.fa")
    script:
        "scripts/make_fasta_file.py"


rule run_beast2:
    input:
        opj(PHYLO_DIR, "mutations.fa")
    output:
        directory("output/beast2"),
        "output/beast2/beast2.model.RDS",
        "output/beast2/tree.log",
    log:
        opj(PHYLO_LOG_DIR, "beast2.log")
    params:
        RSCRIPT=RSCRIPT,
        beast2_dir=BEAST2_DIR,
        site_model=BEAST2_SITE_MODEL,
        clock_model=BEAST2_CLOCK_MODEL,
        tree_prior=BEAST2_TREE_PRIOR,
        iterations=BEAST2_ITERATIONS,
        sample_interval=BEAST2_SAMPLE_INTERVAL,
        rng_seed=BEAST2_RNG_SEED,
    shell:
        """{RSCRIPT} scripts/run_beast2.R -i {input} -o {output[0]} \
            --beast2_dir {params.beast2_dir} \
            --site_model {params.site_model} \
            --clock_model {params.clock_model} \
            --tree_prior {params.tree_prior} \
            --iterations {params.iterations} \
            --sample_interval {params.sample_interval} \
            --rng_seed {params.rng_seed} >& {log}
        """


rule summarize_beast2:
    input:
        "output/beast2/beast2.model.RDS",
    output:
        opj(PHYLO_DIR, "summary.txt"),
        opj(PHYLO_DIR, "summary.ess.txt"),
    log:
        opj(PHYLO_LOG_DIR, "summary.log")
    params:
        RSCRIPT=RSCRIPT,
        sample_interval=BEAST2_SAMPLE_INTERVAL,
        burnin=BEAST2_BURNIN,
    shell:
        """{params.RSCRIPT} scripts/summarize_beast2.R \
            {input[0]} {params.sample_interval} {params.burnin} \
            {output[0]} {output[1]} \
            >& {log}
        """


perc_burnin = int(round(float(BEAST2_BURNIN) / BEAST2_ITERATIONS * 100))
assert perc_burnin > 0 and perc_burnin < 100

rule combine_trees:
    input:
        "output/beast2/tree.log",
    output:
        opj(PHYLO_DIR, "beast2.trees.nexus.txt")
    log:
        opj(PHYLO_LOG_DIR, "beast2.trees.nexus.log")
    params:
        LOGCOMBINER=LOGCOMBINER,
        perc_burnin=perc_burnin,
    shell:
        """
        {params.LOGCOMBINER} -b {params.perc_burnin} \
            -log output/beast2/tree.log -o {output} >& {log}
        """
        

rule make_mcc_tree:
    input:
        opj(PHYLO_DIR, "beast2.trees.nexus.txt")
    output:
        opj(PHYLO_DIR, "max_clade_cred.nexus.txt")
    log:
        opj(PHYLO_LOG_DIR, "max_clade_cred.nexus.log")
    params:
        JAVA=JAVA,
        BEAST2_DIR=BEAST2_DIR,
    shell:
        # -burnin 0 because burning already removed by tree combiner
        """
        {params.JAVA} -Xms1g -Xmx32g \
            -Dlauncher.wait.for.exit=true \
            -Duser.language=en \
            -Djava.library.path={params.BEAST2_DIR}/lib \
            -cp {params.BEAST2_DIR}/lib/launcher.jar \
            beast.app.treeannotator.TreeAnnotatorLauncher \
            -burnin 0 \
            -heights mean \
            -lowMem \
            {input} {output} >& {log}
        """


rule analyze_mcc_tree:
    input:
        opj(PHYLO_DIR, "max_clade_cred.nexus.txt"),
        opj(DEMUX_DIR, "cells.txt"),
    output:
        opj(PHYLO_DIR, "max_clade_cred.newick.txt"),
        opj(PHYLO_DIR, "max_clade_cred.rerooted.newick.txt"),
        opj(PHYLO_DIR, "max_clade_cred.dist.txt"),
        opj(PHYLO_DIR, "max_clade_cred.metadata.txt"),
        opj(PHYLO_DIR, "max_clade_cred.rerooted.metadata.txt"),
    log:
        opj(PHYLO_LOG_DIR, "max_clade_cred.analysis.log")
    params:
        RSCRIPT=RSCRIPT,
    shell:
        """{params.RSCRIPT} scripts/analyze_mcc_tree.R \
            {input[0]} {input[1]} \
            {output[0]} {output[1]} {output[2]} {output[3]} {output[4]}
            >& {log}
        """


def calc_model_indexes(
    min_perc, max_perc, total_iterations, sample_interval):
    # Return a tuple of (i_min_tree, i_max_tree).  i_min_tree and
    # i_max_tree are indexes into model$tree or model$estimates.
    # Indexes are 1-based and inclusive indexes into the <model>$trees
    # variable.
    import math
    
    assert min_perc >= 0 and min_perc <= 100
    assert max_perc > 1 and max_perc <= 100
    assert min_perc < max_perc

    # <model>$trees goes from 1:(num_samples+1) (1-based, inclusive,
    # like R).  trees[1] is sample 0.
    #num_samples = total_iterations / sample_interval
    #i_min = 1
    #i_max = num_samples + 1
    
    # Calculate the ideal starting iteration.
    x1 = total_iterations * min_perc / 100.0
    x2 = total_iterations * max_perc / 100.0
    # Convert to the number of sample intervals.
    x1 = int(math.floor(x1 / sample_interval))
    x2 = int(math.floor(x2 / sample_interval))
    #if x2-x1 < min_samples:
    #    return None
    # Convert to 1-based, inclusive indexes.  Also take into account
    # that the first sample is iteration 0.
    x1 = x1 + 1
    x2 = x2 + 1
    return x1, x2

zoom_min, zoom_max = 90, 100  # how much to zoom in.
x = calc_model_indexes(
    zoom_min, zoom_max, BEAST2_ITERATIONS, BEAST2_SAMPLE_INTERVAL)
i_min, i_max = x
assert i_max-i_min+1 >= 2, "Not enough samples for parameter plots"


rule plot_model_estimates:
    input:
        "output/beast2/beast2.model.RDS"
    output:
        opj(PHYLO_DIR, "posterior.pdf"),
        opj(PHYLO_DIR, "posterior_%03d_%03d.pdf" % (zoom_min, zoom_max)),
        opj(PHYLO_DIR, "tree_height.pdf"),
        opj(PHYLO_DIR, "yule_model.pdf"),      # for yule model only
    log:
        opj(PHYLO_LOG_DIR, "plot_model_estimates.log")
    params:
        RSCRIPT=RSCRIPT,
        BEAST2_BURNIN=BEAST2_BURNIN,
        i_zoom_min=i_min,
        i_zoom_max=i_max,
    shell:
        """{params.RSCRIPT} scripts/plot_model_estimates.R \
            {input[0]} {params.i_zoom_min} {params.i_zoom_max} \
            {output[0]} {output[1]} {output[2]} {output[3]} \
            {params.BEAST2_BURNIN} \
            >& {log}
        """


rule plot_mcc_tree:
    input:
        opj(PHYLO_DIR, "max_clade_cred.newick.txt"),
        opj(DEMUX_DIR, "cells.txt"),
    output:
        opj(PHYLO_DIR, "tree.mcc.pdf"),
    log:
        opj(PHYLO_LOG_DIR, "tree.mcc.log")
    params:
        RSCRIPT=RSCRIPT,
        color_tree_by_category=True,
        tree_height=8,
        tree_width=6,
        tree_layout="rectangular",
    shell:
        """{params.RSCRIPT} scripts/plot_tree.R \
            {input[0]} {input[1]} \
            {params.color_tree_by_category} \
            {params.tree_height} {params.tree_width} \
            {params.tree_layout} {output[0]} \
            >& {log}
        """


perc_min, perc_max = 90, 100  # how much to zoom in.
x = calc_model_indexes(
    zoom_min, zoom_max, BEAST2_ITERATIONS, BEAST2_SAMPLE_INTERVAL)
i_min, i_max = x

MAX_SAMPLES = 100   # Don't plot too many trees, or will be messy.
if i_max-i_min+1 > MAX_SAMPLES:
    i_min = max(i_max-MAX_SAMPLES, 0)
 

rule plot_densitree:
    input:
        "output/beast2/beast2.model.RDS"
    output:
        opj(PHYLO_DIR, "densitree.pdf"),
    log:
        opj(PHYLO_LOG_DIR, "densitree.log")
    params:
        RSCRIPT=RSCRIPT,
        i_min=i_min,
        i_max=i_max,
        tree_height=8,
        tree_width=6,
        tree_layout="rectangular",
    shell:
        """{params.RSCRIPT} scripts/plot_densitree.R \
            {input[0]} \
            {params.i_min} {params.i_max} \
            {params.tree_height} {params.tree_width} \
            {params.tree_layout} {output[0]} \
            >& {log}
        """
