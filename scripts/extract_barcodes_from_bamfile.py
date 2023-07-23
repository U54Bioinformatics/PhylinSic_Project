# Functions:
# hash_var
# find_bam_files
# read_cell_file
# clean_cell_file
# is_all_channel_1
# get_barcodes_from_sample
# write_bam_header
# make_barcode_file
## grep_barcode_from_text_file
## grep_barcode_from_bam_file

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


def find_bam_files(bam_dir):
    import os
    x = os.listdir(bam_dir)
    x = [x for x in x if x.endswith(".bam")]
    x = [os.path.join(bam_dir, x) for x in x]
    return x


def read_cell_file(cell_file, samples):
    # Return a tuple of (<cells>, <used_samples>).  cells is a set of
    # cells in format: <sample>_<barcode>, and may be cleaned up
    # versions of those found in the cell_file.
    #
    # <barcode> includes only the DNA bases, and will include the "-1"
    # GEM chip channel at the end.  If no GEM channel is given, will
    # assume "-1".
    #
    # <sample> will match one of the <used_samples>.  <used_samples>
    # is the subset of <samples> with a cell from the cell file.
    assert samples
    
    # Make a list of the samples from the CellRanger files.
    # User's samples may be hashed versions of these.
    sample2i = {}  # (possibly hashed) sample -> index into samples
    for i, sample in enumerate(samples):
        h1 = sample
        h2 = hash_var(sample)
        h3 = hash_var(sample, can_start_with_number=True)
        sample2i[h1] = i
        sample2i[h2] = i
        sample2i[h3] = i

    # Read in each of the cells and clean up.
    x = open(cell_file).read()
    x = x.split()
    x = [x.strip() for x in x]
    x = [x for x in x if x]
    cells = set(x)

    # Check to make sure the format of good_cells look reasonable.
    DNA_BASES = set("ACGT")
    good_cells = set()
    used_samples = set()
    for cell in cells:
        # If user accidentally added ".bam" to name of cell, remove
        # it.
        if cell.endswith(".bam"):
            cell = cell[:-4]
        x = cell.rsplit("_", 1)
        assert len(x) == 2, \
            "Cell not formatted as <sample>_<barcode>: %s" % cell
        sample, barcode = x
        clean_barcode = barcode
        # Parse out the GEM chip channel.
        x = barcode.split("-")
        assert len(x) <= 2
        clean_barcode = x[0]
        gem_channel = "1"
        if len(x) >= 2:
            gem_channel = x[1]
        assert set(clean_barcode).issubset(DNA_BASES), \
            "Barcode contains non-DNA bases in cell: %s" % cell

        x = ["%s_<barcode>" % x for x in sorted(sample2i)]
        x1 = "\n".join(x[:5])
        x2 = "%s_%s" % (samples[0], barcode)
        x = "%s_<barcode>" % sample
        assert sample in sample2i, (
            'Sample from cell file ("%s") doesn\'t match samples from '
            "CellRanger output.\n%s\nCell names should look like: %s" % (
                x, x1, x2))
        clean_sample = samples[sample2i[sample]]
        clean_cell = "%s_%s-%s" % (clean_sample, clean_barcode, gem_channel)
        good_cells.add(clean_cell)
        used_samples.add(clean_sample)
    return good_cells, used_samples


def clean_cell_file(infile, bam_files, outfile):
    import os
    
    # Make sure sample names match.
    # <sample>.bam
    x = [os.path.split(x)[1] for x in bam_files]
    x = [os.path.splitext(x)[0] for x in x]
    sample_names = x

    # data/cellranger.out/006_R_LUNG/outs/possorted_genome_bam.bam
    good_cells, used_samples = read_cell_file(infile, sample_names)
    assert good_cells, "No cells found: %s" % demux_file
    assert used_samples, "No samples found: %s" % demux_file
    # Write out the clean file.
    x = ["%s\n" % x for x in sorted(good_cells)]
    open(outfile, 'w').writelines(x)

    return good_cells, used_samples


def is_all_channel_1(cells):
    for x in cells:
        x = x.rsplit("_", 1)
        assert len(x) == 2
        x = x[1].split("-")
        assert len(x) == 2
        if x[1] != "1":
            return False
    return True


def get_barcodes_from_sample(cells, sample):
    # cells include cells across any sample.  Will return a set of the
    # barcodes for the cells from this sample.
    
    barcodes = set()
    for cell in cells:
        x = cell.strip().rsplit("_", 1)
        assert len(x) == 2, "Cell not <sample>_<barcode> format: %s" % cell
        # cells should be clean and verified.        
        s, b = x
        if s != sample:
            continue
        assert b.find("-") >= 0  # make sure has GEM well
        barcodes.add(b)
    return barcodes


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


def make_barcode_file(barcodes, outfile):
    # According to 10x, grep on CB:Z (for cellranger) or BX:Z (for
    # longranger).
    # https://kb.10xgenomics.com/hc/en-us/articles/360022448251-Is-there-way-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
    handle = open(outfile, 'w')
    for b in barcodes:
        # CB:Z:GCCAAATTCACATACG-1
        assert b.find("-") >= 0  # make sure includes GEM well
        x1 = "CB:Z:%s" % b
        x2 = "BX:Z:%s" % b
        print(x1, file=handle)
        print(x2, file=handle)
    handle.close()


#def grep_barcode_from_text_file(text_file, barcode_file, outfile, log_file):
#    import subprocess
#    import pipes
#    
#    x1 = "cat %s" % pipes.quote(text_file)
#    x2 = "LC_ALL=C grep -F -f %s" % pipes.quote(barcode_file)
#    x3 = pipes.quote(outfile)
#    x4 = pipes.quote(log_file)
#    x = "%s | %s 1> %s 2> %s" % (x1, x2, x3, x4)
#    subprocess.run(x, close_fds=True, check=True, shell=True)
#
#
#def grep_barcode_from_bam_file(bam_file, barcode_file, outfile, log_file):
#    import subprocess
#    import pipes
#    
#    # samtools view $BAM_FILE | LC_ALL=C grep -F -f barcode_file > outfile
#    x1 = "samtools view %s" % pipes.quote(bam_file)
#    x2 = "LC_ALL=C grep -F -f %s" % pipes.quote(barcode_file)
#    x3 = pipes.quote(outfile)
#    x4 = pipes.quote(log_file)
#    x = "%s | %s 1> %s 2> %s" % (x1, x2, x3, x4)
#    subprocess.run(x, close_fds=True, check=True, shell=True)


def main():
    import os
    
    cellranger_dir = snakemake.input[0]
    cell_file_messy = snakemake.input[1]
    out_dir = snakemake.output[0]
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    temp_dir = "temp/split_cells_by_batch"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    cell_file = os.path.join(temp_dir, "cell.txt")

    bam_files = find_bam_files(cellranger_dir)
    # Make sure cells match what is in the CellRanger BAM files.
    # e.g. BRST002_003_INF_ANT_CHEST_MASS_ATTCTTGTCTTGTGCC-1
    x = clean_cell_file(cell_file_messy, bam_files, cell_file)
    good_cells, used_samples = x
    # If all the GEM chip channels are -1, then strip them.
    strip_channel = is_all_channel_1(good_cells)


    # Use 16 cells per batch.  When demultiplexing, need to open
    # this number of file handles at once.  So don't make it too
    # big or there will be too many file handles.
    BATCH_SIZE = 16
    jobs = []
    for orig_bam_file in bam_files:
        x = os.path.split(orig_bam_file)[1]
        sample = os.path.splitext(x)[0]
        if sample not in used_samples:
            continue
        
        header_file = os.path.join(temp_dir, "%s.header" % sample)
        sam_file = os.path.join(temp_dir, "%s.sam" % sample)
        log_file = os.path.join(temp_dir, "%s.log" % sample)
        barcode_file = os.path.join(temp_dir, "%s.barcodes.txt" % sample)

        # Make this barcode file now.  Then split the barcodes
        # into batches.
        x = get_barcodes_from_sample(good_cells, sample)
        barcodes = sorted(x)
 
        batch_num = 0
        for i in range(0, len(barcodes), BATCH_SIZE):
            batch_num += 1
            batch_barcodes = barcodes[i:i+BATCH_SIZE]
            root = "%s.%03d" % (sample, batch_num)
            batch_barcode_file = os.path.join(
                out_dir, "%s.barcodes.txt" % root)
            batch_sam_file = os.path.join(temp_dir, "%s.sam" % root)
            batch_log_file = os.path.join(temp_dir, "%s.log" % root)
            demux_script_file = os.path.join(temp_dir, "%s.py" % root)
            demux_log_file = os.path.join(temp_dir, "%s.log" % root)

            x = type("", (), {})()
            x.sample = sample
            x.orig_bam_file = orig_bam_file
            x.barcodes = barcodes
            x.barcode_file = barcode_file
            x.sam_file = sam_file
            x.log_file = log_file
            x.header_file = header_file
            x.batch_barcodes = batch_barcodes
            x.batch_barcode_file = batch_barcode_file
            x.batch_sam_file = batch_sam_file
            x.batch_log_file = batch_log_file
            x.demux_script_file = demux_script_file
            x.demux_log_file = demux_log_file
            jobs.append(x)
        assert jobs

    # Make the headers for each samples.
    seen = set()
    for j in jobs:
        if j.header_file in seen:
            continue
        seen.add(j.header_file)
        write_bam_header(j.orig_bam_file, j.header_file)

    # Make the barcode file for each sample.
    seen = set()
    for j in jobs:
        if j.barcode_file in seen:
            continue
        make_barcode_file(j.barcodes, j.barcode_file)

    # Samples are split into batches.
    # Make the barcode file for each batch.
    for j in jobs:
        make_barcode_file(j.batch_barcodes, j.batch_barcode_file)

#    # As an optimization, generate smaller versions of the
#    # <bam_files> first that contain only the data from the cells
#    # of interest.
#    seen = set()
#    for j in jobs:
#        if j.sam_file in seen:
#            continue
#        seen.add(j.sam_file)
#        grep_barcode_from_bam_file(
#            j.orig_bam_file, j.barcode_file, j.sam_file, j.log_file)
# 
##         # Now pull out the reads for each batch of barcodes.
##         all_commands = []
##         commands = []
##         for j in jobs:
##             x = _grep_barcode_file_cmd(
##                 j.sam_file, False, j.batch_barcode_file, j.batch_sam_file,
##                 j.batch_log_file)
##             all_commands.append(x)
##             if not filelib.exists_nz(j.batch_sam_file):
##                 commands.append(x)
##         parallel.pshell(commands, num_slots=num_cores, hours=1)
## 
##         # Now demultiplex with the python script.
## 
##         # Previous implementation sometimes generated bad SAM files.
##         # Got errors with "samtools flagstats" saying that file is
##         # truncated.  (samtools quickcheck is fine).  Doesn't always
##         # happen, but when it does, seems to affect many sam files, as
##         # if all jobs were killed prematurely.  Not sure why?
##         # Addressed by removing buffering, and also used smaller sets
##         # of barcodes so don't need to open and close files.
## 
##         #x = mlib.calc_nodes_and_cores(
##         #    len(sample_dirs), num_cores, max_threads=8)
##         #nc, nt = x
## 
##         # Demultiplex each of the alignments.
##         hours = []
##         commands = []
##         for j in jobs:
##             hdr_map = singularitylib.Bind(j.header_file)
##             sam_map = singularitylib.Bind(j.batch_sam_file)
##             out_map = singularitylib.Bind(out_dir)
##             x = singularitylib.run_cmd(
##                 "python2", "python", hdr_map, sam_map, out_map)
##             cmd = singularitylib.quote_and_join(x)
##             params = {
##                 "SAMPLE" : j.sample,
##                 "BARCODES" : j.batch_barcodes,
##                 "STRIP_CHANNEL" : strip_channel,
##                 "HEADER_FILE" : hdr_map.container_path,
##                 "ALIGN_FILE" : sam_map.container_path,
##                 "OUT_DIR" : out_map.container_path,
##                 }
##             x = scriptlib.make_py_cmd(
##                 cmd, DEMULTIPLEX_PY, j.demux_script_file, j.demux_log_file,
##                 params)
##             commands.append(x)
##             # 30 Gb BAM file takes ~5 hours.  Not sure where the
##             # 96 hours came from.  Seems like overkill.
##             # SAM file should be shorter, since it's uncompressed.
##             h = 96
##             GB = 1024 * 1024 * 1024
##             size = filelib.filesize(j.batch_sam_file)
##             if size < 32 * GB:
##                 h = 8
##             elif size < 4 * GB:
##                 h = 1
##             hours.append(h)
##         memory = 8
##         parallel.pshell(
##             commands, num_slots=num_cores, hours=hours, memory=memory)
##         metadata["num_cores"] = num_cores
## 
##         # Make sure all the scripts completed successfully.
##         x = [j.demux_log_file for j in jobs]
##         filelib.assert_exists_nz_many(x)
## 
##         # Read the number of cells found.
##         results = []
##         for j in jobs:
##             x = open(j.demux_log_file).read().strip()
##             num_cells = int(x)
##             results.append(num_cells)
## 
##         cells_found = sum(results)
##         msg = "No cells demultiplexed."
##         msg = "%s  Maybe samples don't match demux_cell_file." % msg
##         assert cells_found, msg
    
    
if __name__ == '__main__':
    main()



## from Module import AbstractModule
## 
## class Module(AbstractModule):
##     def __init__(self):
##         AbstractModule.__init__(self)
## 
##     def run(
##         self, network, antecedents, out_attributes, user_options, num_cores,
##         out_dir):
##         import os
##         from genomicode import parallel
##         from genomicode import filelib
##         from genomicode import filefinder
##         from genomicode import scriptlib
##         from genomicode import singularitylib
##         from Betsy import module_utils as mlib
## 
##         bam_node = antecedents
##         bam_filenames = filefinder.find_bam_files(
##             bam_node.identifier, find_cram=True)
##         assert bam_filenames, "No .bam files."
##         filelib.safe_mkdir(out_dir)
##         metadata = {}
## 
##         # <bam_node.identifier>/
##         #   <sample>.bam                                <orig_bam_file>
##         # <out_dir>/
##         #   <sample>_<barcode>.sam
##         # <temp_dir>/
##         #   demux_file.txt         <clean_demux_file>
##         #   <sample>.header        <header_file>   # SAM header for this sample
##         #   <sample>.barcodes.txt  <barcode_file>  # barcodes of good cells
##         #   <sample>.sam           <sam_file>      # only cells of interest
##         #   <sample>.log           <log_file>
##         #   <sample>.<batch>.barcodes.txt  <batch_barcode_file>
##         #   <sample>.<batch>.sam           <batch_sam_file>
##         #   <sample>.<batch>.log           <batch_log_file>
##         #   <sample>.<batch>.py            <script_file>
##         #   <sample>.<batch>.log           <log_file>
##         opj = os.path.join
##         temp_dir = "temp"
##         filelib.safe_mkdir(temp_dir)
##         clean_demux_file = opj(temp_dir, "demux_file.txt")
## 
##         # Algorithm:
##         # 1.  Make a <sam_file> from <orig_bam_file> that contains
##         #     only the cells of interest.
##         #     grep for reads matching <barcode_file>.
##         # 2.  Split the <sam_file> into batches of 16 cells each.
##         # 3.  Use a Python script to deconvolute each batch into
##         #     individual cells.
##         
##         # SAM files for single cells range from 20 Mb (40k reads) to 6
##         # Gb (10m reads).
##         
##         # Can't just match the barcode because multiple samples may
##         # have the same barcode.
## 
##         # cells should be in the format:
##         # <sample>_<barcode>  (may include -1 GEM well).
##         #
##         # Don't try to be too clever and match partial sample names.
##         # There are problems with this, e.g.
##         # User's sample looks like:                007_<barcode>
##         # Cell Ranger sample looks like:   BRST007_005_<barcode>
## 
##         demux_file = mlib.get_user_option(
##             user_options, "demux_cell_file", check_file=True, not_empty=True)
## 
##         # Clean up the demultiplex file.  Make sure sample names
##         # match.
##         sample_names = [mlib.splitpath(x)[1] for x in bam_filenames]
##         good_cells, used_samples = _read_demux_file(
##             demux_file, sample_names)
##         assert good_cells, "No cells found: %s" % demux_file
##         assert used_samples, "No samples found: %s" % demux_file
##         # Write out the clean file.
##         x = ["%s\n" % x for x in sorted(good_cells)]
##         open(clean_demux_file, 'w').writelines(x)
## 
##         # cells in good_cells look like:
##         # BRST002_003_INF_ANT_CHEST_MASS_ATTCTTGTCTTGTGCC-1
## 
##         # If all the GEM chip channels are -1, then strip them.
##         strip_channel = True
##         for x in good_cells:
##             x = x.rsplit("_", 1)
##             assert len(x) == 2
##             x = x[1].split("-")
##             assert len(x) == 2
##             if x[1] != "1":
##                 strip_channel = False
##         
##         # Use 16 cells per batch.  When demultiplexing, need to open
##         # this number of file handles at once.  So don't make it too
##         # big, or there will be too many file handles.
##         BATCH_SIZE = 16
##         jobs = []
## 
##         for orig_bam_file in bam_filenames:
##             x, sample, x = mlib.splitpath(orig_bam_file)
##             if sample not in used_samples:
##                 continue
##             header_file = opj(temp_dir, "%s.header" % sample)
##             sam_file = opj(temp_dir, "%s.sam" % sample)
##             log_file = opj(temp_dir, "%s.log" % sample)
##             barcode_file = opj(temp_dir, "%s.barcodes.txt" % sample)
## 
##             # Make this barcode file now.  Then split the barcodes
##             # into batches.
##             x = _get_barcodes_from_sample(good_cells, sample)
##             barcodes = sorted(x)
## 
##             batch_num = 0
##             for i in range(0, len(barcodes), BATCH_SIZE):
##                 batch_num += 1
##                 batch_barcodes = barcodes[i:i+BATCH_SIZE]
##                 root = "%s.%03d" % (sample, batch_num)
##                 batch_barcode_file = opj(temp_dir, "%s.barcodes.txt" % root)
##                 batch_sam_file = opj(temp_dir, "%s.sam" % root)
##                 batch_log_file = opj(temp_dir, "%s.log" % root)
##                 demux_script_file = opj(temp_dir, "%s.py" % root)
##                 demux_log_file = opj(temp_dir, "%s.log" % root)
##                 
##                 x = filelib.GenericObject(
##                     sample=sample,
##                     orig_bam_file=orig_bam_file,
##                     barcodes=barcodes,
##                     barcode_file=barcode_file,
##                     sam_file=sam_file,
##                     log_file=log_file,
##                     header_file=header_file,
##                     batch_barcodes=batch_barcodes,
##                     batch_barcode_file=batch_barcode_file,
##                     batch_sam_file=batch_sam_file,
##                     batch_log_file=batch_log_file,
##                     demux_script_file=demux_script_file,
##                     demux_log_file=demux_log_file,
##                     )
##                 jobs.append(x)
##         assert jobs
## 
##         # Make the headers for each sample.
##         commands = []
##         seen = set()
##         for j in jobs:
##             if j.header_file in seen:
##                 continue
##             seen.add(j.header_file)
##             if not filelib.exists_nz(j.header_file):
##                 x = _write_bam_header, (j.orig_bam_file, j.header_file), {}
##                 commands.append(x)
##         parallel.pyfun(commands, num_slots=num_cores)
## 
##         # Make the barcode files for each sample and batch.
##         commands = []
##         seen = set()
##         for j in jobs:
##             if j.barcode_file not in seen:
##                 seen.add(j.barcode_file)
##                 if not filelib.exists_nz(j.barcode_file):
##                     x = _make_barcode_file, (j.barcodes, j.barcode_file), {}
##                     commands.append(x)
##             if not filelib.exists_nz(j.batch_barcode_file):
##                 x = j.batch_barcodes, j.batch_barcode_file
##                 x = _make_barcode_file, x, {}
##                 commands.append(x)
##         parallel.pyfun(commands, num_slots=num_cores)
## 
##         # As an optimization, generate smaller versions of the
##         # <bam_files> first that contain only the data from the cells
##         # of interest.
##         all_commands = []
##         commands = []
##         seen = set()
##         for j in jobs:
##             if j.sam_file in seen:
##                 continue
##             seen.add(j.sam_file)
##             x = _grep_barcode_file_cmd(
##                 j.orig_bam_file, True, j.barcode_file, j.sam_file, j.log_file)
##             all_commands.append(x)
##             if not filelib.exists_nz(j.sam_file):
##                 commands.append(x)
##         parallel.pshell(commands, num_slots=num_cores, hours=4)
## 
##         # Now pull out the reads for each batch of barcodes.
##         all_commands = []
##         commands = []
##         for j in jobs:
##             x = _grep_barcode_file_cmd(
##                 j.sam_file, False, j.batch_barcode_file, j.batch_sam_file,
##                 j.batch_log_file)
##             all_commands.append(x)
##             if not filelib.exists_nz(j.batch_sam_file):
##                 commands.append(x)
##         parallel.pshell(commands, num_slots=num_cores, hours=1)
## 
##         # Now demultiplex with the python script.
## 
##         # Previous implementation sometimes generated bad SAM files.
##         # Got errors with "samtools flagstats" saying that file is
##         # truncated.  (samtools quickcheck is fine).  Doesn't always
##         # happen, but when it does, seems to affect many sam files, as
##         # if all jobs were killed prematurely.  Not sure why?
##         # Addressed by removing buffering, and also used smaller sets
##         # of barcodes so don't need to open and close files.
## 
##         #x = mlib.calc_nodes_and_cores(
##         #    len(sample_dirs), num_cores, max_threads=8)
##         #nc, nt = x
## 
##         # Demultiplex each of the alignments.
##         hours = []
##         commands = []
##         for j in jobs:
##             hdr_map = singularitylib.Bind(j.header_file)
##             sam_map = singularitylib.Bind(j.batch_sam_file)
##             out_map = singularitylib.Bind(out_dir)
##             x = singularitylib.run_cmd(
##                 "python2", "python", hdr_map, sam_map, out_map)
##             cmd = singularitylib.quote_and_join(x)
##             params = {
##                 "SAMPLE" : j.sample,
##                 "BARCODES" : j.batch_barcodes,
##                 "STRIP_CHANNEL" : strip_channel,
##                 "HEADER_FILE" : hdr_map.container_path,
##                 "ALIGN_FILE" : sam_map.container_path,
##                 "OUT_DIR" : out_map.container_path,
##                 }
##             x = scriptlib.make_py_cmd(
##                 cmd, DEMULTIPLEX_PY, j.demux_script_file, j.demux_log_file,
##                 params)
##             commands.append(x)
##             # 30 Gb BAM file takes ~5 hours.  Not sure where the
##             # 96 hours came from.  Seems like overkill.
##             # SAM file should be shorter, since it's uncompressed.
##             h = 96
##             GB = 1024 * 1024 * 1024
##             size = filelib.filesize(j.batch_sam_file)
##             if size < 32 * GB:
##                 h = 8
##             elif size < 4 * GB:
##                 h = 1
##             hours.append(h)
##         memory = 8
##         parallel.pshell(
##             commands, num_slots=num_cores, hours=hours, memory=memory)
##         metadata["num_cores"] = num_cores
## 
##         # Make sure all the scripts completed successfully.
##         x = [j.demux_log_file for j in jobs]
##         filelib.assert_exists_nz_many(x)
## 
##         # Read the number of cells found.
##         results = []
##         for j in jobs:
##             x = open(j.demux_log_file).read().strip()
##             num_cells = int(x)
##             results.append(num_cells)
## 
##         cells_found = sum(results)
##         msg = "No cells demultiplexed."
##         msg = "%s  Maybe samples don't match demux_cell_file." % msg
##         assert cells_found, msg
## 
##         return metadata
## 
## 
##     def name_outfile(self, antecedents, out_attributes, user_options):
##         return "sam"
## 
## 
## def _get_barcodes_from_sample(cells, sample):
##     # cells include cells across any sample.  Will return a set of the
##     # barcodes for the cells from this sample.
##     
##     barcodes = set()
##     for cell in cells:
##         x = cell.strip().rsplit("_", 1)
##         assert len(x) == 2, "Cell not <sample>_<barcode> format: %s" % cell
##         # cells should be clean and verified.        
##         s, b = x
##         if s != sample:
##             continue
##         assert b.find("-") >= 0  # make sure has GEM well
##         barcodes.add(b)
##     return barcodes
##     
## 
## 
## def _make_barcode_file(barcodes, outfile):
##     # According to 10x, grep on CB:Z (for cellranger) or BX:Z (for
##     # longranger).
##     # https://kb.10xgenomics.com/hc/en-us/articles/360022448251-Is-there-way-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
##     handle = open(outfile, 'w')
##     for b in barcodes:
##         # CB:Z:GCCAAATTCACATACG-1
##         assert b.find("-") >= 0  # make sure includes GEM well
##         x1 = "CB:Z:%s" % b
##         x2 = "BX:Z:%s" % b
##         print >>handle, x1
##         print >>handle, x2
##     handle.close()
## 
## 
## def _grep_barcode_file_cmd(
##         bam_file, is_sam_or_bam, barcode_file, outfile, log_file):
##     # Return a string.
##     # If not is_sam_or_bam, then just a text file to grep.
##     from genomicode import singularitylib
##     from genomicode import filelib
## 
##     # samtools view $BAM_FILE | LC_ALL=C grep -F -f barcode_file > outfile
##     if is_sam_or_bam:
##         F = singularitylib.File
##         x = singularitylib.run_cmd("samtools", "samtools", "view", F(bam_file))
##         x1 = singularitylib.quote_and_join(x)
##     else:
##         x1 = "cat %s" % filelib.quote(bam_file)
##     x2 = "LC_ALL=C grep -F -f %s" % filelib.quote(barcode_file)
##     x3 = filelib.quote(outfile)
##     x4 = filelib.quote(log_file)
##     x = "%s | %s 1> %s 2> %s" % (x1, x2, x3, x4)
##     return x
##     
## 
## 
## DEMULTIPLEX_PY = """
## 
## def main():
##     import os
##     from genomicode import samtools
## 
##     MAX_HANDLES = 32
##     # If there's 8 batches and 32 handles, can open files for 256
##     # cells.  There may be 1000-5000 cells for each 10x run.  Chances
##     # of cache miss are pretty high, if cells accessed randomly.
##     # In practice, cache misses are pretty rare (< 1:1000).
##     SAMPLE = {{SAMPLE}}
##     BARCODES = {{BARCODES}}
##     STRIP_CHANNEL = {{STRIP_CHANNEL}}  # Strip the channel from the filename.
##     HEADER_FILE = {{HEADER_FILE}}
##     ALIGN_FILE = {{ALIGN_FILE}}
##     OUT_DIR = {{OUT_DIR}}
## 
##     barcodes = set(BARCODES)
##     # Don't want too many open files.
##     assert len(barcodes) <= 64, "Too many barcodes.  Not enough batches."
## 
##     # Open a file for each barcode.
##     barcode2handle = {}
##     for barcode in barcodes:
##         x = barcode
##         if STRIP_CHANNEL:
##             assert x[-2] == "-"
##             x = x[:-2]
##         filename = os.path.join(OUT_DIR, "%s_%s.sam" % (SAMPLE, x))
##         barcode2handle[barcode] = open(filename, 'w', buffering=0)
## 
##     # Write the header to each of the outfiles.
##     header = list(open(HEADER_FILE))
##     for handle in barcode2handle.itervalues():
##         handle.writelines(header)
## 
##     found = set()
##     for line in open(ALIGN_FILE):
##         align = samtools.parse_alignment(line)
##         
##         # CR  Cell barcode
##         # CY  Cell barcode read quality
##         # CB  Cell barcode, error-corrected, confirmed.  
##         #     <barcode>-<gem_chip_channel>
##         # UR  UMI
##         # UY  UMI read quality
##         # UB  UMI, error-corrected
## 
##         # According to 10x, may be BX:Z for longranger.
##         # CB:Z:CCACGAGCATCAACCA-1
##         x = align.tags.get("CB", None)
##         if not x:   # no barcode
##             # Ignore reads without a barcode.
##             continue
##         x, barcode = x
##         assert x == "Z"
##         #if barcode.endswith("-1"):
##         #    barcode = barcode[:-2]
##         # This should not happen if we had grep'd properly.
##         if barcode not in barcode2handle:
##             continue
##         barcode2handle[barcode].write(line)
##         found.add(barcode)
## 
##     # Print out the number of cells (not number of reads).
##     print len(found)
## 
## main()
## 
## """
## 
## 
