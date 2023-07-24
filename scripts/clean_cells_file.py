# Functions:
# hash_var
# read_cell_file

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


def read_cell_file(cell_file, samples):
    # Return a tuple of:
    #   (<cells>, <cell2category>, <cell2outgroup>, <used_samples>).
    # 
    # cells is a set of cells in format: <sample>_<barcode>, and may
    # be cleaned up versions of those found in the cell_file.
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
    cell2category = {}
    cell2outgroup = {}
    i_cell = None
    i_category = None
    i_outgroup = None
    cells = set()
    for line in open(cell_file):
        cols = line.rstrip("\r\n").split("\t")
        if i_cell is None:
            assert "Cell" in cols, \
                'File is missing a column with header "Cell".'
            i_cell = cols.index("Cell")
            if "Category" in cols:
                i_category = cols.index("Category")
            if "Outgroup" in cols:
                i_outgroup = cols.index("Outgroup")
            continue
        cell = cols[i_cell].strip()
        if not cell:
            continue
        cells.add(cell)
        if i_category is not None:
            cell2category[cell] = cols[i_category]
        if i_outgroup is not None:
            cell2outgroup[cell] = cols[i_outgroup]
    
    # Check to make sure the format of good_cells look reasonable.
    DNA_BASES = set("ACGT")
    good_cells = set()
    good_cell2category = {}
    good_cell2outgroup = {}
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
        if cell in cell2category:
            good_cell2category[clean_cell] = cell2category[cell]
        if cell in cell2outgroup:
            good_cell2outgroup[clean_cell] = cell2outgroup[cell]
    return good_cells, good_cell2category, good_cell2outgroup, used_samples


def main():
    import os
    
    user_cell_file = snakemake.input[0]
    bam_files = snakemake.input[1:]
    clean_cell_file = snakemake.output[0]
    
    #temp_dir = "temp/demux_cellranger"
    #if not os.path.exists(temp_dir):
    #    os.makedirs(temp_dir)
    #cell_file = os.path.join(temp_dir, "cell.txt")

    # Make sure cells match what is in the CellRanger BAM files.  They
    # should look like:
    #   BRST002_003_INF_ANT_CHEST_MASS_ATTCTTGTCTTGTGCC-1

    # Parse out the sample names from the BAM files.
    #   <sample>.bam
    x = [os.path.split(x)[1] for x in bam_files]
    x = [os.path.splitext(x)[0] for x in x]
    sample_names = x

    x = read_cell_file(user_cell_file, sample_names)
    good_cells, cell2category, cell2outgroup, used_samples = x
    assert good_cells, "No cells found: %s" % user_cell_file
    assert used_samples, "No samples found: %s" % user_cell_file
    
    # Write out the clean file.
    handle = open(clean_cell_file, 'w')
    header = "Cell", "Category", "Outgroup"
    print("\t".join(header), file=handle)
    for cell in sorted(good_cells):
        category = cell2category.get(cell, "")
        outgroup = cell2outgroup.get(cell, "")
        x = cell, category, outgroup
        assert len(x) == len(header)
        print("\t".join(x), file=handle)
    handle.close()

    
if __name__ == '__main__':
    main()
