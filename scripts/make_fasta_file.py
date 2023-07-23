def choose_het(ref, alt):
    POSSIBLE = ["A", "C", "G", "T"]
    for x in POSSIBLE:
        if x not in [ref, alt]:
            return x
    raise AssertionError


def write_fasta(title, sequence, width=60, handle=None):
    # title does not include ">".
    import sys

    handle = handle or sys.stdout
    w = handle.write
    w(f">{title}\n")
    i = 0
    while i < len(sequence):
        s = sequence[i:i+width]
        i += width
        w(f"{s}\n")


def main():
    genotype_file = snakemake.input[0]
    out_file = snakemake.output[0]

    # Write out the bases for each cell as FASTA.  Don't change
    # the order of the bases.

    # Read everything into memory.  Hopefully, it all fits.
    inhandle = open(genotype_file, encoding="utf-8")
    x = inhandle.readline()
    header = x.rstrip("\r\n").split("\t")
    cell_names = header[1:]

    # Read as a list of strings.  (Should take less memory than a
    # list of lists.)
    bases = []
    for x in inhandle:
        cols = x.rstrip("\r\n").split("\t")
        assert len(cols) == len(header)

        # Convert the genotypes to bases.
        site = cols[0]
        x = site.split("_")
        assert len(x) == 4
        chrom, pos, ref, alt = x
        assert len(ref) == 1, x
        assert len(alt) == 1, x
        geno2base = {
            "" : "-",   # IUPAC uses a single dash or hyphen for gaps.
            "R" : ref,
            "A" : alt,
            "H" : choose_het(ref, alt),
            }
        x = cols[1:]
        x = [geno2base.get(x, None) for x in x]
        assert None not in x
        x = "".join(x)
        assert len(x) == len(cell_names)
        bases.append(x)
    assert bases, "No data"
    inhandle.close()

    handle = open(out_file, 'w', encoding="utf-8")
    for i in range(len(cell_names)):
        x = [x[i] for x in bases]
        seq = "".join(x)
        write_fasta(cell_names[i], seq, handle=handle)
    handle.close()


if __name__ == '__main__':
    main()
