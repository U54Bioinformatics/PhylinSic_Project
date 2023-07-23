def main():
    import os
    import extract_sample_barcodes as barcodelib
    
    cell_file = snakemake.input[0]
    barcode_file = snakemake.output[0]
    sample = snakemake.params.sample
    # is a string, because used to construct file names.
    batch_num = int(snakemake.params.batch)
    batch_size = snakemake.params.batch_size

    # Read in a set of the cells.
    cells = barcodelib.read_cell_file(cell_file)
    
    # Make a list of the barcodes that are seen in sample.
    x = barcodelib.get_barcodes_from_sample(cells, sample)
    barcodes = sorted(x)   # barcodes should have GEM well.

    # Find the barcodes associated with this batch.
    i = batch_num * batch_size
    batch_barcodes = barcodes[i:i+batch_size]
    
    barcodelib.make_barcode_file(batch_barcodes, barcode_file)
    
    
if __name__ == '__main__':
    main()
