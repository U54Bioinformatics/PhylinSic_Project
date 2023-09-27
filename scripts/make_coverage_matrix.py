def iter_coord(filename):
    # Yield lines that have the same Chrom, Pos, Ref, Alt

    # File should have headers:
    # Chrom  Pos  Ref  Alt  Sample  Coverage  Cell

    it = open(filename)
    x = next(it)
    header = x.rstrip("\r\n").split("\t")
    assert len(header) == 7, f"Invalid header in file: {filename}"
    assert header == [
        "Chrom", "Pos", "Ref", "Alt", "Sample", "Coverage", "Cell"]

    lines = []
    prev_coord = None
    for line in it:
        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) == 7, f"Invalid line in file: {filename}\n{line}"
        coord = cols[:4]
        if coord != prev_coord:
            if lines:
                yield lines
            lines = []
            prev_coord = coord
        lines.append(cols)
    if lines:
        yield lines
            
        
def main():
    import os
    
    in_file = snakemake.input[0]
    pb_matrix_file = snakemake.input[1]
    out_file = snakemake.output[0]

    # in_file:
    # Chrom  Pos  Ref  Alt  Sample  Coverage  Cell

    # Use the Ref and Alt alleles from the pseudobulk.
    handle = open(pb_matrix_file, encoding="utf-8")
    x0 = handle.readline()
    x1 = handle.readline()
    x2 = handle.readline()
    header2 = x2.rstrip("\r\n").split("\t")
    assert header2[:4] == ["Chrom", "Pos", "Ref", "Alt"]
    coord2ra = {}   # (chrom, pos (int)) -> (ref, alt)
    for x in handle:
        cols = x.rstrip("\r\n").split("\t")
        chrom, pos, ref, alt = cols[:4]
        pos = int(pos)
        coord2ra[(chrom, pos)] = (ref, alt)
    handle.close()

    # Make two passes through in_file.  First, make a list of all the
    # cells.  Then, pull out the coverage information.
    cells = set()
    for lines in iter_coord(in_file):
        for cols in lines:
            assert len(cols) == 7
            cells.add(cols[6])
    cells = sorted(cells)

    outhandle = open(out_file, 'w', encoding="utf-8")
    header = ["Chrom", "Pos", "Ref", "Alt"] + cells
    print("\t".join(header), file=outhandle)

    for lines in iter_coord(in_file):
        #chrom, pos, ref, alt = lines[0][:4]
        chrom, pos = lines[0][:2]
        pos = int(pos)
        assert (chrom, pos) in coord2ra
        ref, alt = coord2ra[(chrom, pos)]
        
        cov_cols = [""] * len(cells)
        for cols in lines:
            coverage = cols[5]
            cell = cols[6]
            i_cell = cells.index(cell)
            cov_cols[i_cell] = coverage
        x = [chrom, pos, ref, alt] + cov_cols
        assert len(x) == len(header)
        print("\t".join(map(str, x)), file=outhandle)

    outhandle.close()


if __name__ == '__main__':
    main()
