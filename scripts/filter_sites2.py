def main():
    import bisect

    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    num_bases = int(snakemake.params.remove_close_sites)
    assert num_bases is None or num_bases >= 0

    inhandle = open(in_file, encoding="utf-8")
    outhandle = open(out_file, 'w', encoding="utf-8")

    x = inhandle.readline()
    header = x.rstrip("\r\n").split("\t")
    assert header
    # Chrom  Pos  Ref  Alt  <site>  ...
    assert len(header) > 4
    assert header[:4] == ["Chrom", "Pos", "Ref", "Alt"]
    print("\t".join(header), file=outhandle)

    coords = []  # list of (chrom, pos (int)), sorted
    for line in inhandle:
        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) == len(header)

        chrom, pos, ref, alt = cols[:4]
        pos = int(pos)

        # Arbitrarily keep the first one seen in the file.
        i = bisect.bisect_left(coords, (chrom, pos))
        i1 = max(i-2, 0)
        i2 = min(i+2, len(coords))
        is_close = False
        for j in range(i1, i2):
            c, p = coords[j]
            if c != chrom:
                continue
            if abs(pos-p) > num_bases:
                continue
            is_close = True
            break
        coords.insert(i, (chrom, pos))

        if is_close:
            continue
        print(line, file=outhandle, end="")

    inhandle.close()
    outhandle.close()


if __name__ == '__main__':
    main()
