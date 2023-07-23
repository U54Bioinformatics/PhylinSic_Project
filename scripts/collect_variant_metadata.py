def make_metadata_file(variant_file, outfile):
    # Make a first pass to get an overview of the information
    # inside.  Should make this fast, since file is big.

    # Format:
    # Coord   <chrom>  <pos>  <ref>  <alt>
    # Source  <caller>  <source>
    # Sample  <sample>

    coord_data = set()
    samples = set()
    caller2source = set()  # (caller, source)

    # Want Caller, Sample, and Source.
    outhandle = open(outfile, 'w', encoding="utf-8")

    handle = open(variant_file, encoding="utf-8")
    x = next(handle)
    assert x, "Empty file: %s" % variant_file
    header = x.rstrip("\r\n").split("\t")
    i_chrom = header.index("Chrom")
    i_pos = header.index("Pos")
    i_ref = header.index("Ref")
    i_alt = header.index("Alt")
    i_caller = header.index("Caller")
    i_sample = header.index("Sample")
    i_source = header.index("Source")
    for line in handle:
        if line is None:
            break
        x = line.rstrip("\r\n").split("\t")
        assert len(x) == len(header)
        # Get the coordinate.
        chrom, pos, ref, alt = x[i_chrom], x[i_pos], x[i_ref], x[i_alt]
        y = chrom, int(pos), ref, alt
        if y not in coord_data:
            coord_data.add(y)
            y = ["Coord"] + list(y)
            print("\t".join(map(str, y)), file=outhandle)
        # Get the other information.
        caller = x[i_caller]
        sample = x[i_sample]
        source = x[i_source]
        if (caller, source) not in caller2source:
            caller2source.add((caller, source))
            y = "Source", caller, source
            print("\t".join(y), file=outhandle)
        if sample not in samples:
            samples.add(sample)
            y = "Sample", sample
            print("\t".join(y), file=outhandle)
    handle.close()


def main():
    import os

    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    make_metadata_file(in_file, out_file)


if __name__ == '__main__':
    main()

