# Functions:
# read_cell_file
# get_barcodes_from_sample
# make_barcode_file


def read_cell_file(filename):
    # Return a set of cells.  The file should be cleaned up.  Cells in
    # format: <sample>_<barcode>-<gem_well>
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
    return cells


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


def main():
    import os
    
    cell_file = snakemake.input[0]
    barcode_file = snakemake.output[0]
    sample = snakemake.params.sample

    # Read in a set of the cells.
    cells = read_cell_file(cell_file)

    # Make a list of the barcodes that are seen in sample.
    x = get_barcodes_from_sample(cells, sample)
    barcodes = sorted(x)   # barcodes should have GEM well.
    
    # Make the barcode file.
    make_barcode_file(barcodes, barcode_file)
    
    
if __name__ == '__main__':
    main()
