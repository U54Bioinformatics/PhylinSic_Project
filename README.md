# PhylinSic

This contains code for phylogenetic inference from single-cell RNA-Seq
data.  This is described in:

Liu X, Griffiths J, Bishara I, Liu J, Bild AH, and Chang JT.
Phylogenetic inference from single-cell RNA-seq data.
bioRxiv 2022.09.27.509725.  doi:
https://doi.org/10.1101/2022.09.27.509725


# Installation

This pipeline depends on software that must be pre-installed on your
system.  The dependencies are listed below.  I have also included
instructions for installing with conda.  However, everyone's
environment is a little different and software changes over time, so
you may need to do things differently.

* [snakemake](https://snakemake.readthedocs.io/en/stable/)

* samtools

```
conda install -c bioconda samtools
```

* picard

```
conda install -c bioconda picard
```


* GATK4

```
conda install -c bioconda gatk4
```

* Varscan

```
conda install -c bioconda varscan
```

* numpy

```
conda install -c anaconda numpy
```


* scipy

```
conda install -c anaconda scipy
```

* R/babette

```
conda install -c bioconda r-argparse
conda install -c bioconda bioconductor-ggtree
conda install -c conda-forge r-phangorn
conda install -c conda-forge r-treetools

Rscript -e 'install.packages("babette", repos="http://cran.us.r-project.org")
```


* BEAST2

```
conda install -c bioconda beast2
```


Remember when I mentioned that you may need to do things differently?
Well, I needed to.  I ran into an incompatibility between Java and R
that I could not resolve.  The rJava installed by conda did not
work.

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

This is a problem because R and Java are needed to run BEAST2, the
last step of the pipeline.  I ended up installing R+libraries, Java,
and BEAST2 outside of Conda, and then configuring the last step of the
pipeline to use my manually installed version.  If you run into
similar problems, there are instructions in the Snakefile for how to
force the last rule to use a specific installation of R and BEAST2.

Another issue: the picard and gatk4 packages work with different
versions of openjdk.  To work aroud this, I use picard from conda and
configure `gatk` to use an external installation.


# Running the pipeline

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

Snakefile is downloaded from this repository.

`cells.txt` is a text file where each line gives the name of a cell to
be analyzed.  This should include only the high-quality cells in the
data set (e.g. preprocess the data with CellRanger, filter for cells
based on read depth, mitochondria, droplet analysis, etc.)  The names
of the cells should be in the format <sample>_<barcode>, where the
<sample> should match the name of a sample analyze by CellRanger.  An
example file can be found [here](cells.txt).

The `cellranger` directory contains the output from the CellRanger
preprocessing.  Each subdirectory of `cellranger` should contain the
data for a sample.  Although CellRanger creates many files, the only
ones we are about are the BAM files `possorted_genome_bam.bam` that
contain the alignments for samples.  All the other files are ignored.
Please do not move the BAM files, or we won't be able to find them.

`genome.fa` is the FASTA-formatted file containing the reference
genome that the BAM files were aligned to.  This needs to be the same
genome that CellRanger used.

`known_sites1.vcf.gz` should contain the mutation sites XXX


2.  Configure Snakefile.

Open up Snakefile in your favorite text editor and configure the
parameters inside.  This is how you control how the inference is done.
The parameters are described in the file.  They are set to reasonable
defaults, although you may need to customize them to your data set
(e.g. if you are analyzing tumors with a very large number of
mutations and need stricter mutation filtering.)


3.  Run snakemake.

While in <your directory>, start snakemake by running something like:

```
snakemake --cores 8 --latency-wait 0 --rerun-triggers mtime
```

Hopefully, this will run all the way through without errors.  If you
run into bugs in the pipeline, please let me know.  XXX


4.  Get results from XXX.

```
<your directory>/
    output/
        beast2/
```







# Pipeline

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
XXX get rid of duplicate SNPs
keep_only_snps                   We only want SNPs, remove INDELs.
parse_vcf_files                  Pull out variant information from VCF files.
collect_variant_metadata         Collect XXX metadata about variants.
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
merge_coverage_by_batch          Merge each XXX batch.
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



# Directory structure:

data/
  cells.txt       XXX
  cellranger/     A directory with all the results from CellRanger.
                  It should be organized according to the following
                  structure:
                  cellranger_results/
                    <sample1>/
                      <sample>.mri
                      _cmdline
                      _filelist
                      outs/
                        possorted_genome_bam.bam
                        possorted_genome_bam.bam.bai
                        ...
                      SC_RNA_COUNTER_CS/
                        ...
                    <sample2>/
                      ... as in <sample1>
                    <sample3>/
                      ... as in <sample1>
  genome.fa       Reference genome used by CellRanger for alignment.
  known_sites1.vcf.gz
  known_sites1.vcf.gz.tbi
  known_sites2.vcf.gz
  known_sites2.vcf.gz.tbi
  known_sites3.vcf.gz
  known_sites3.vcf.gz.tbi
output/
  # 1.  Prepare the reference files.
  genome.fa
  genome.fa.fai
  genome.dict
  # 2.  Demultiplex CellRanger output into single cells.
  cells.txt             # A cleaned up version of the cells file.
  02_demux/{sample}.bam
  02_demux/{sample}.barcodes.txt
  02_demux/{sample}.header.txt
  02_demux/{sample}.alignments.txt
  02_demux/{sample}.{batch}.barcodes.txt
  02_demux/{sample}.{batch}.alignments.txt
  02_demux/{sample}.{batch}.cells.sam/
  # 3.  Preprocess the single cell files.
  03_preproc/{sample}.{batch}.cells.bam/
  03_preproc/{sample}.{batch}.cells.rg.bam/
  03_preproc/{sample}.{batch}.cells.rg_contig.bam/
  03_preproc/{sample}.{batch}.cells.rg_contig_index.bam/
  03_preproc/{sample}.{batch}.cells.rg_contig_index_snt.bam/
  03_preproc/{sample}.{batch}.cells.rg_contig_index_snt.report/
  03_preproc/{sample}.{batch}.cells.rg_contig_index_snt_recal.bam/
  # 4.  Merge single cells to a pseudobulk sample.
  04_pbulk/{sample}.{batch}.merged.bam
  04_pbulk/{sample}.{batch}.bam
  04_pbulk/{sample}.merged.bam
  04_pbulk/{sample}.bam
  04_pbulk/pseudobulk.merged.bam
  04_pbulk/pseudobulk.cleaned.bam
  04_pbulk/pseudobulk.bam
  04_pbulk/pseudobulk.bam.bai
  # 5.  Call variants on the pseudobulk sample.
  05_call/interval_{interval}.intervals
  05_call/variants.{interval}.vcf
  05_call/variants.snp_only.{interval}.vcf
  05_call/variants.table.txt
  05_call/variants.metadata.txt
  05_call/variants.unfiltered.matrix.{batch}.txt
  05_call/variants.unfiltered.matrix.txt
  05_call/variants.matrix.txt
  05_call/variants.coord.txt
  # 6.  Make a site x cell matrix containing the ref/alt read depth.
  06_matrix/{sample}.{batch}.raw.pileup/
  06_matrix/{sample}.{batch}.pileup/
  06_matrix/{sample}.{batch}.coverage.vcf/
  06_matrix/{sample}.{batch}.coverage.txt/
  06_matrix/{sample}.{batch}.coverage.with_cells.txt
  06_matrix/{sample}.{batch}.coverage.merged.txt
  06_matrix/coverage.table.txt
  06_matrix/coverage.matrix.txt
  # 7.  Filter the sites.
  07_filter/coverage.matrix.{batch}.txt
  07_filter/coverage.matrix.filter1.{batch}.txt
  07_filter/coverage.matrix.filter1.txt
  07_filter/coverage.matrix.filter2.txt
  07_filter/coverage.metadata.txt
  07_filter/ref_count.all.txt
  07_filter/alt_count.all.txt
  07_filter/knn.features.txt
  07_filter/ref_count.features.txt
  07_filter/alt_count.features.txt
  07_filter/knn.neighbors.{batch}.txt
  07_filter/knn.neighbors.txt
  07_filter/ref_count.{batch}.txt
  07_filter/alt_count.{batch}.txt
  07_filter/genotypes.{batch}.txt
  07_filter/probabilities.{batch}.txt
  07_filter/genotypes.txt
  07_filter/probabilities.txt
  07_filter/genotype.count.txt
  07_filter/mixed_genotypes.txt
  07_filter/coverage.matrix.filter3.txt
  07_filter/coverage.matrix.filter4.txt
  # 8.  Call/impute genotypes.
  08_genotype/coverage.metadata.txt
  08_genotype/ref_count.all.txt
  08_genotype/ref_count.all.txt
  08_genotype/knn.features.txt
  08_genotype/ref_count.features.txt
  08_genotype/alt_count.features.txt
  08_genotype/knn.neighbors.0.txt
  08_genotype/knn.neighbors.txt
  08_genotype/genotypes.txt
  08_genotype/probabilities.txt
  # 9.  Make the phylogeny.
  mutations.fa
  beast2/
    beast2.infile.txt
    beast2.outfile.txt
    beast2.model.RDS
    screen.log
    trace.log
    tree.log
logs/
  genome.fa.fai.log
  genome.dict.log
  {sample}.alignments.log
  {sample}.{batch}.alignments.log
  {sample}.{batch}.cells.rg.log
  {sample}.{batch}.cells.rg_contig.log
  {sample}.{batch}.cells.rg_contig_index.log
  {sample}.{batch}.cells.rg_contig_index_snt.log
  {sample}.{batch}.cells.rg_contig_index_snt_report.log
  {sample}.{batch}.cells.rg_contig_index_snt_recal.log
  # 4.  Merge single cells to a pseudobulk sample.
  {sample}.{batch}.pb.log
  {sample}.pb.log
  {sample}.pb.clean.log
  pseudobulk.merged.log
  pseudobulk.cleaned.log
  pseudobulk.sort.log
  pseudobulk.index.log
  # 5.  Call variants on the pseudobulk sample.
  interval_{interval}.intervals
  pseudobulk.clean.contig.{interval}.log
temp/
XXX RENAME THESE
  {sample}.pb.files
  {sample}.pb.merged.header
  {sample}.pb.header
  {sample}.{batch}.pb.header
  {sample}.{batch}.pb.clean.header
  {sample}.{batch}.pb.files
  pseudobulk.pb.files
  pseudobulk.merged.header
  pseudobulk.header
  {sample}.{batch}.coverage.merged.01.txt
  {sample}.{batch}.coverage.merged.02.txt
  {sample}.{batch}.coverage.merged.03.txt
  coverage.01.txt
  coverage.02.txt
  coverage.03.txt
scripts/
  demux_cellranger.py
  XXX




# FAQ

How much disk space?

How long it takes to run?



