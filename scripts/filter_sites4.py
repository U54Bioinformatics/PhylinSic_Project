def main():
    mutation_file = snakemake.input[0]
    out_file = snakemake.output[0]

    from filter_sites1 import parse_rav

    x = snakemake.params.keep_n_sites_seen_in_most_cells
    assert x is None or x > 0
    keep_n_sites = x

    outhandle = open(out_file, 'w', encoding="utf-8")


    # Figure out the number of sites to keep if
    # keep_n_sites_seen_in_most_cells was requested.
    ns_cutoff = None
    # Get the distribution of the number of sites each cell has.
    nsites2count = {}
    if keep_n_sites is not None:
        inhandle = open(mutation_file, encoding="utf-8")
        x = inhandle.readline()
        header = x.rstrip("\r\n").split("\t")
        assert header
        # Chrom  Pos  Ref  Alt  <site>  ...
        assert len(header) > 4
        assert header[:4] == ["Chrom", "Pos", "Ref", "Alt"]

        for x in inhandle:
            cols = x.rstrip("\r\n").split("\t")
            assert len(cols) == len(header)

            ns = 0
            for x in cols[4:]:
                if not x:
                    continue
                ref, alt, vaf = parse_rav(x)
                if ref > 0 or alt > 0:
                    ns += 1
            nsites2count[ns] = nsites2count.get(ns, 0) + 1
        assert nsites2count, "No sites"

        # Figure out the nsites cutoff.
        # Sort by decreasing number of sites.
        schwartz = [(-ns, ns, count) for (ns, count) in nsites2count.items()]
        schwartz.sort()
        counts = [x[1:] for x in schwartz]

        n = 0
        while n < keep_n_sites and counts:
            ns, count = counts.pop(0)
            n += count
            ns_cutoff = ns
        assert ns_cutoff is not None

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

        if ns_cutoff is not None:
            ns = 0
            for x in cols[4:]:
                if not x:
                    continue
                ref, alt, vaf = parse_rav(x)
                if ref > 0 or alt > 0:
                    ns += 1
            if ns < ns_cutoff:
                continue

        print(line, file=outhandle, end="")

    inhandle.close()
    outhandle.close()


if __name__ == '__main__':
    main()
