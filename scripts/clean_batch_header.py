# Merging BAM files with samtools can result in duplicated lines in
# the header.  This can be a problem if the header exceeds the maximum
# limit for a BAM file (2Gb, 2^31).

# @HD        1
# @SQ       84
# @RG      255    One line for each cell that makes up this file.
# @PG   12,240    Programs run on each cell.
# @CO  118,730    Merged from each cell.  Lots of duplicates.
# Getting rid of duplicates will get it down to around 10% of its size.

def dedup_header(in_file, header_file, clean_header_file, out_file, log_file):
    import subprocess
    import extract_bam_header

    # Write out a cleaned up version of the header.
    extract_bam_header.write_bam_header(in_file, header_file)
    handle = open(clean_header_file, 'w')
    seen = set()
    for line in open(header_file):
        if line in seen:
            continue
        seen.add(line)
        handle.write(line)
    handle.close()

    cmd = [
        "samtools", "reheader", "-P", clean_header_file, in_file,
        ]
    out_handle = open(out_file, 'w')
    log_handle = open(log_file, 'w')
    subprocess.run(
        cmd, stdout=out_handle, stderr=log_handle, close_fds=True,
        check=True)
    out_handle.close()
    log_handle.close()



def main():
    in_file = snakemake.input[0]
    out_file = snakemake.output[0]
    header_file = snakemake.params.header_file
    clean_header_file = snakemake.params.clean_header_file
    log_file = snakemake.params.log_file

    dedup_header(
        in_file, header_file, clean_header_file, out_file, log_file)



    

    
    
if __name__ == '__main__':
    main()
