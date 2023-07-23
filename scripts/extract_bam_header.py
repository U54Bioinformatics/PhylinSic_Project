# Functions:
# write_bam_header


def write_bam_header(bam_file, outfile):
    import subprocess

    # samtools view bam_file -H >& outfile
    cmd = [
        "samtools",
        "view",
        bam_file,
        "-H",
        ]
    outhandle = open(outfile, 'w', encoding="utf-8")
    subprocess.run(
        cmd, stdout=outhandle, stderr=subprocess.STDOUT, close_fds=True,
        check=True)
    outhandle.close()


def main():
    import os
    
    bam_file = snakemake.input[0]
    out_file = snakemake.output[0]

    write_bam_header(bam_file, header_file)
    
    
if __name__ == '__main__':
    main()
