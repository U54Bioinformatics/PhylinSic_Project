
def main():
    import os

    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    # Write a positions file.
    # <chrom> <pos 1-based>

    # GATK uses 1-based coordinates.  Keep in that way.

    outhandle = open(out_file, 'w')
    handle = open(in_file)
    header1 = next(handle)  # Read the header
    header2 = next(handle)
    header3 = next(handle)
    cols = header3.rstrip("\r\n").split("\t")
    assert len(cols) == 6
    assert cols[:2] == ["Chrom", "Pos"]
    seen = set()
    for line in handle:
        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) == 6
        chrom, pos = cols[:2]
        pos = int(pos)
        if (chrom, pos) in seen:
            continue
        seen.add((chrom, pos))
        x = [chrom, pos]
        print("\t".join(map(str, x)), file=outhandle)

    handle.close()
    outhandle.close()
        


if __name__ == '__main__':
    main()

