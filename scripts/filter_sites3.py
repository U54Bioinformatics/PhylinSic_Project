def main():
    mutation_file = snakemake.input[0]
    good_sites_file = snakemake.input[1]
    out_file = snakemake.output[0]

    # Read in the good_sites_file.
    positions = set()
    handle = open(good_sites_file, 'r', encoding="utf-8")
    x = handle.readline()
    header = x.rstrip("\r\n").split("\t")
    assert len(header) == 4
    assert header == ["Chrom", "Pos", "Ref", "Alt"]
    for x in handle:
        cols = x.rstrip("\r\n").split("\t")
        assert len(cols) == 4
        positions.add(tuple(cols))
    handle.close()


    outhandle = open(out_file, 'w', encoding="utf-8")

    # Read and filter the file.
    inhandle = open(mutation_file, encoding="utf-8")
    x = inhandle.readline()
    header = x.rstrip("\r\n").split("\t")
    assert header
    # Chrom  Pos  Ref  Alt  <site>  ...
    assert len(header) > 4
    assert header[:4] == ["Chrom", "Pos", "Ref", "Alt"]
    print("\t".join(header), file=outhandle)

    for line in inhandle:
        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) == len(header)

        x = tuple(cols[:4])
        if x not in positions:
            continue

        print(line, file=outhandle, end="")

    inhandle.close()
    outhandle.close()


if __name__ == '__main__':
    main()
