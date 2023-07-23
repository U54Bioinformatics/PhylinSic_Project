
def main():
    import os

    in_file = snakemake.input[0]
    out_file = snakemake.output[0]
    
    min_total_reads = snakemake.params.keep_calls_with_min_total_reads
    min_alt_reads = snakemake.params.keep_calls_with_min_alt_reads
    min_vaf = snakemake.params.keep_calls_with_min_vaf

    if min_total_reads is not None:
        assert min_total_reads >= 0 and min_total_reads < 1E6
    if min_alt_reads is not None:
        assert min_alt_reads >= 0 and min_alt_reads < 10000
    if min_vaf is not None:
        assert abs(min_vaf-1.0) > 1E-4, \
               "Minimum allele frequency cannot be 1."
        assert min_vaf >= 0.0 and min_vaf < 1.0, \
               "Allele frequency should be between 0 and 1.  min vaf is %s" % \
               min_vaf
    
    assert not (
        min_total_reads is None and
        min_alt_reads is None and
        min_vaf is None)

    outhandle = open(out_file, 'w', encoding="utf-8")

    #                           Num Callers  <Sample>
    #                                        HaplotypeCaller (RNA)
    # Chrom   Pos   Ref   Alt   <Sample>     Ref/Alt/VAF
    handle = open(in_file)
    header1 = next(handle)
    header2 = next(handle)
    header3 = next(handle)
    assert header1
    assert header2
    assert header3
    header1 = header1.rstrip("\r\n").split("\t")
    header2 = header2.rstrip("\r\n").split("\t")
    header3 = header3.rstrip("\r\n").split("\t")
    assert len(header1) == 6
    assert len(header2) == 6
    assert len(header3) == 6

    print("\t".join(header1), file=outhandle)
    print("\t".join(header2), file=outhandle)
    print("\t".join(header3), file=outhandle)

    for line in handle:
        cols = line.rstrip("\r\n").split("\t")
        x = cols[-1].split("/")
        assert len(x) == 3
        num_ref, num_alt, vaf = int(x[0]), int(x[1]), float(x[2])
        if min_alt_reads is not None and num_alt < min_alt_reads:
            continue
        if min_total_reads is not None and num_alt + num_ref < min_total_reads:
            continue
        if min_vaf is not None and vaf < min_vaf:
            continue
        print(line, end="", file=outhandle)

    handle.close()
    outhandle.close()


if __name__ == '__main__':
    main()

