# Algorithm.
# 1.  Pull out the headers for each of the BAM files.
#     Make sure that samples don't conflict.
# 2.  Merge into batches of MAX_MERGE_FILES.
#     samtools will give an error if we try to merge too many files at
#     once.
# 3.  Remove duplicates in headers.
#     When merging samples, the concatenated headers might have
#     duplicated lines.
#     If merging many small files (e.g. scRNA-Seq), may end up with a
#     header that exceeds the BAM header limit of 2*31 bytes.
# 4.  Merge the batches into the final file.
# 5.  Remove duplicates in headers.


# Functions:
# hash_var
# read_cell_file
# clean_cell_file

def hash_var(name, can_start_with_number=False):
    import re
    # Fix the header to be a python variable.
    x = str(name)
    # Replace all non-word character with _.
    x = re.sub(r"\W", "_", x)
    if not can_start_with_number:
        # Replace initial numbers with Xnumber.
        x = re.sub(r"^(\d)", r"X\1", x)
    return x


def main():
    import os
    
    bam_files = snakemake.input
    pseudobulk_file = snakemake.output[0]



    
    
if __name__ == '__main__':
    main()
