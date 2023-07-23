def main():
    import os
    
    cellranger_dir = snakemake.input[0]
    out_file = snakemake.output[0]
    sample = snakemake.params.sample

    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # <cellranger_dir>/<sample>/outs/possorted_genome_bam.bam
    x = os.path.join(
        cellranger_dir, sample, "outs", "possorted_genome_bam.bam")
    filename = os.path.realpath(x)
    assert os.path.exists(filename), \
           f"I could not find BAM file for sample {sample}: {filename}"
    os.symlink(filename, out_file)
    
    
if __name__ == '__main__':
    main()
