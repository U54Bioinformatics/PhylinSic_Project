This repository contains code for phylogenetic inference from
single-cell RNA-Seq data that is described in:

Liu X, Griffiths J, Bishara I, Liu J, Bild AH, and Chang JT.
Phylogenetic inference from single-cell RNA-seq data.
bioRxiv 2022.09.27.509725.  doi:
https://doi.org/10.1101/2022.09.27.509725

The pipeline is implemented with Snakemake and takes 10X Genomics
Chromium single cell RNA-Seq data processed with CellRanger and
generates a phylogenetic tree.

If you would like to do your own processing, we also provide just the
R functions that will smooth, impute, and call the genotypes for
scRNA-Seq data as described in the paper.  You can get that code by
downloading the [smoothing.R](smoothing.R) file.

This README contains the following sections:
- [Installation](#Installation)
- [Quick Start](#Start)
- [Pipeline](#Pipeline)
- [FAQ](#FAQ)


# <A NAME="Installation">Installation</A>

This pipeline requires the following software to be installed on your
system.  I have included instructions for installing these with conda.
However, everyone's environment is a little different, and software
changes over time, so you may need to do things differently to get
things to work.

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)

* Samtools
* Picard
* GATK4
* Varscan
* Python3 + numpy + scipy
* R + babette + other libraries
* BEAST2


```
conda install -c bioconda samtools
conda install -c bioconda picard
conda install -c bioconda gatk4
conda install -c bioconda varscan

conda install -c anaconda numpy
conda install -c anaconda scipy

conda install -c bioconda r-argparse
conda install -c r r-rcolorbrewer
conda install -c conda-forge r-ggplot2
conda install -c bioconda bioconductor-ggtree
conda install -c conda-forge r-phangorn
conda install -c conda-forge r-treetools

Rscript -e 'install.packages("babette", repos="http://cran.us.r-project.org")

conda install -c bioconda beast2
```

*Note*: Using Conda, I ran into an incompatibility between Java and R
that I could not resolve.  First, the rJava installed by Conda was not
linked properly.

```
> library(rJava)
Error: package or namespace load failed for ‘rJava’:
 .onLoad failed in loadNamespace() for 'rJava', details:
  call: dyn.load(file, DLLpath = DLLpath, ...)
  error: unable to load shared object '.../R/library/rJava/libs/rJava.so':
  libjvm.so: cannot open shared object file: No such file or directory
```

After I fixed this issue (either by setting LD_LIBRARY_PATH or
encoding a path with -rpath), I still could not run rJava:

```
> rJava::.jinit()
* runs forever *
```

To work around this problem, I installed R+libraries, Java, and BEAST2
outside of Conda, and then configured the pipeline to use my version
instead.  If you run into similar problems, there are instructions in
the Snakefile that describe how to do this.

*Another issue*: the picard and gatk4 packages in Conda require
different versions of openjdk.  When I installed gatk4, conda
downgraded openjdk, which broke picard.  To work around this, I used
picard from conda and installed Java+GATK4 manually (and configured
the Snakefile accordingly).


# <A NAME="Start">Quick Start</A>

1.  Set up the input files.

You should create a local directory and set it up like:

```
<your directory>/
    Snakefile
    data/
        cells.txt
        cellranger/
            <sample1>/
                outs/
                    possorted_genome_bam.bam
                    possorted_genome_bam.bam.bai
                    ...
                ...
            <sample2>/
              ... as in <sample1>
            <sample3>/
              ... as in <sample1>
        genome.fa
        known_sites1.vcf.gz
        known_sites1.vcf.gz.tbi
        known_sites2.vcf.gz
        known_sites2.vcf.gz.tbi
        known_sites3.vcf.gz
        known_sites3.vcf.gz.tbi
```

`Snakefile` is downloaded from this repository.

`cells.txt` is a tab-delimited text file.  The first row should
contain headers "Cell", "Outgroup", and "Category".  Each line gives
the name of a cell to be analyzed.  This should include only the
high-quality cells in the data set (e.g. preprocess the data with
CellRanger, filter for cells based on read depth, mitochondria,
droplet analysis, etc.)  The names of the cells should be in the
format \<sample\>_\<barcode\>, where the \<sample\> should match the
name of a sample analyze by CellRanger.  The "Outgroup" and "Category"
columns are optional (you may leave them out of the file), but can
provide extra information that is used when interpreting the
phylogenies.  "Outgroup" is used to indicate (with "yes" or "no")
which cells comprise the outgroup (e.g. the non-cancer cells in a
tumor).  "Category" can be filled with some grouping of cells
(e.g. the cell type, which tumor the cells came from, etc.) and is
only used for plotting.  An example file can be found
[here](cells.txt).

The `cellranger` directory contains the output from the CellRanger
preprocessing.  Each subdirectory of `cellranger` should contain the
data for a sample.  Although CellRanger creates many files, the only
ones we are about are the BAM files `possorted_genome_bam.bam` that
contain the alignments for samples.  All the other files are ignored.
Please do not move the BAM files, or we won't be able to find them.

`genome.fa` is the FASTA-formatted file containing the reference
genome that the BAM files were aligned to.  This needs to be the same
genome that CellRanger used.

`known_sites1.vcf.gz` (and `known_sites2.vcf.gz` and
`known_sites3.vcf.gz`) should contain known mutation sites to be used
for realignment and base quality recalibration.

We typically use:
```
Mills_and_1000G_gold_standard.indels.b37.vcf.gz
1000G_phase1.indels.b37.vcf.gz
dbsnp_138.b37.vcf.gz
```

You can get these files from the Broad online.  Make sure they are
compressed and indexed with bgzip and tabix.




2.  Configure Snakefile.

Download [Snakefile](Snakefile) to the computer you will run the
pipeline on.  Open a text editor and configure the parameters in the
file.  The parameters are set to reasonable defaults, although you may
need to customize them to your data set (e.g. if you are analyzing
tumors with a very large number of mutations and need stricter
mutation filtering.)


3.  Run snakemake.

While in `<your directory>`, start snakemake by running something
like:

```
snakemake --cores 8 --latency-wait 0 --rerun-triggers mtime
```

Hopefully, this will run all the way through without errors.  If you
run into bugs in the pipeline, please file an
[issue](https://github.com/U54Bioinformatics/PhylinSic_Project/issues).


4.  Get the results.

This pipeline will generate a phylogenetic tree of the cells.  The raw
output from the BEAST2 modeling is available in the `output/beast2/`
directory, and some processed results are available in `output/phylogeny/`.

```
<your directory>/
    output/
        beast2/
            beast2.infile.txt
            beast2.outfile.txt
            beast2.model.RDS
            screen.log
            trace.log
            tree.log
        phylogeny/
            mutations.fa              Alignments used for modeling.

            summary.txt               Summary statistics.
            summary.ess.txt           Effective sample size estimates.
            posterior.pdf             Plot of posterior probabilities.
            posterior_090_100.pdf     Posterior of last 10% of samples.
            tree_height.pdf           Plot of tree height estimates.
            yule_model.pdf            Only if "yule" tree prior was selected.

            beast2.trees.nexus.txt         Estimated trees across sampling.
            max_clade_cred.nexus.txt       Max clade credibility tree, NEXUS.
            max_clade_cred.newick.txt      Max clade credibility tree, Newick.
            max_clade_cred.dist.txt        Pairwise distances.
            max_clade_cred.metadata.txt    Metadata describing tree.
            max_clade_cred.rerooted.newick.txt     Rerooted, if Outgroup given.
            max_clade_cred.rerooted.metadata.txt

            tree.mcc.pdf              Dendrogram of max clade credibility tree.
            densitree.pdf             Dendrogram showing sampled trees.

```




# <A NAME="Pipeline">Pipeline</A>

The pipeline is broken down into nine tasks.

*1.  Prepare the reference genome.*

Do basic pre-processing and indexing of the reference genome needed
for subsequent steps.

```
Rules:
copy_ref_genome                  Copy the reference genome to a local path.
index_ref_genome                 Index the reference genome.
create_ref_genome_dict           Create a sequence dictionary.
```


*2.  Demultiplex CellRanger output into single cells.*

Each sample is split into batches of single cells that are processed
in parallel.  The files for each batch of cells are stored in a
directory.

```
Rules:
extract_cellranger_bam           Find the BAM files in the CellRanger output.
extract_bam_header               Pull out the header from a bam file.
clean_cells_file                 Make a clean version of the user's cells.
extract_sample_barcodes          Pull out the barcodes for a sample.
extract_sample_alignments        Pull out the alignments for a sample.
extract_batch_barcodes           Pull out the barcodes for a sample+batch
extract_batch_alignments         Pull out the alignments for a sample+batch.
demux_one_batch                  Demultiplex alignments into single cells.
```


*3.  Preprocess the single cell files.*

These rules work on directories (batches) of single cell
files.

```
Rules:
convert_sam_to_bam               Generate BAM files for processing.
add_read_groups
sort_bam
index_bam
split_n_trim
make_base_recalibration_report
recalibrate_base_quality_score
```

*4.  Merge single cells to a pseudobulk sample.*

The BAM files for each batch of cells is merged into one BAM file.
The BAM files for each sample is merged into a BAM file.  All BAM
files is merged into a pseudobulk BAM file.  The merging is done
piecemeal so that we do not run into problems when trying to merge too
many files at once.  After each merge, we need to remove duplicate
lines from the BAM headers, otherwise they may get too big (>2^31
bytes).

```
Rules:
merge_cells_to_batch             Merge cells from each batch into a BAM file
clean_batch_header               Clean up the headers after merging.
merge_batches_to_sample          Merge batches into one sample file.
clean_sample_header              Clean up the headers after merging.
merge_samples_to_pseudobulk      Merge samples into one pseudobulk file.
clean_pseudobulk_header          Clean up the headers after merging.
sort_pseudobulk_by_contig        Sort the pseudobulk BAM file by contig.
index_pseudobulk_bam             Index the pseudobulk BAM file.
```


*5.  Call variants on the pseudobulk sample.*

Call variants with GATK and keep the ones that meet a read depth
cutoff.  Make a table of the genomic coordinates with variants.  These
will be the candidate sites to use for the phylogeny reconstruction.

```
Rules:
make_interval_file               Break the genome down into intervals.
call_variants                    Call variants on a single interval.
keep_only_snps                   We only want SNPs, remove INDELs.
parse_vcf_files                  Pull out variant information from VCF files.
collect_variant_metadata         Get metadata (coords, cells) about variants.
convert_to_matrix_by_batch       Reformat variant information to a matrix.
merge_matrix_batches             Merge into a single matrix with Ref/Alt/VAF.
filter_variant_calls             Filter variant calls based on read depth.
extract_variant_coordinates      Pull out the genomic coords of the variants.
```


*6.  Make a site x cell matrix containing the ref/alt read depth.*

At the variant sites found in the pseudobulk, count the
reference and alternate read depth for each of the cells.
Represent as a site x cell matrix of ref/alt/vaf data.

```
Rules:
make_pileup_from_cells           Generate pileups for each cell at each site.
clean_pileup_file                Clean up the pileups.
convert_pileup_to_vcf            Convert the pileup into VCF format.
extract_coverage_from_vcf        Pull out the coverage.
add_cell_names_to_coverage       Add the cell names to the coverage files.
merge_coverage_by_batch          Merge coverage for each batch of cells.
merge_coverage_batches           Merge all coverage into one file.
make_coverage_matrix             Reformat into a site x cell matrix.
```


*7.  Filter the sites.*

There will be too many sites to process, so we select the ones that
are most likely to be informative.  Since different filters require
different types of information (e.g. can be done one site at a time,
or requires information from multiple sites; requires predicted
genotype call; etc.), we do the filtering in stages.

Stage 1: Apply the filters that evaluate each site independently of
other sites.  Since the sites don't depend on each other, these
filters are easy to run in parallel.

Stage 2: Apply the filters that require information from multiple
sites.

Stage 3: Apply the filters that require the genotypes to be called.
We call the genotypes before running these filters.

Stage 4: Apply the filters that should be run last.

```
Rules:
split_matrix_for_filter1         Split coverage matrix to filter in parallel.
filter_sites1                    Do stage 1 filtering.
merge_filter1_matrix             Re-merge the filtered files.
filter_sites2                    Do stage 2 filtering.

prepare_geno_calling_files1      Collect info needed to call genotypes.
select_geno_calling_features1    Select features for genotype calling.
extract_geno_calling_features1   Pull out the features.
calc_neighbor_scores1            Score similarity among all the sites.
select_neighbors1                Select neighbors to use for imputation.
split_ref_matrix1                Split the matrices to call genotypes.
split_alt_matrix1                Split the matrices to call genotypes.
call_genotypes1                  Call or impute the genotypes, in batches.
merge_genotype_files1            Merge the files containing genotype calls.
merge_probabilities_files1       Merge the files with the prob of genotypes.

count_genotypes                  Count the genotypes seen at each site.
list_sites_with_mixed_genotypes  List the sites with a mix of genotypes.
filter_sites3                    Do stage 3 filtering.
filter_sites4                    Do stage 4 filtering.
```

*8.  Call/impute genotypes.*

After the filtering is done, then make the final set of genotype calls
from the (high quality) filtered data.  Because these files have
already been filtered and contain a small number of sites, we do not
need to call the genotypes in batches as above.  We can call on the
whole file, simplying the pipeline.

```
Rules:
prepare_geno_calling_files2
select_geno_calling_features2
extract_geno_calling_features2
calc_neighbor_scores2
select_neighbors2
call_genotypes2
```


*9.  Make the phylogeny.*

Use the called/imputed genotypes and estimate a phylogeny.

```
Rules:
make_fasta_file                  Make a FASTA file with the genotypes.
run_beast2                       Run BEAST2 to generate the phylogenies.
```







# <A NAME="FAQ">FAQ</A>

- What kind of computer do I need to run this pipeline?

In short, you need a lot of everything.  More CPUs help pre-process
the cells faster, faster CPUs helps the phylogenies to be calculated
faster, RAM is needed to load large matrices of mutations into memory,
and disk is needed to process the alignments for all the cells.

How much of everything you need depends on the size of your data set
and how long you are willing to wait for the computation to finish.
But if I had to set a minimum lower bound, let's say you'll need a
LINUX server with 16 cores, 32 Gb RAM, and 1 Tb hard drive.  But I
would recommend 92 cores, 64 Gb RAM, and 4 Tb hard drive.
