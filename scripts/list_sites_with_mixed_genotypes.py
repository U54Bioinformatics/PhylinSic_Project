def list_good_sites(genotype_count_file, perc_mixed_calls, outfile):
    assert perc_mixed_calls > 1 and perc_mixed_calls < 50

    outhandle = open(outfile, 'w', encoding="utf-8")
    outheader = "Chrom", "Pos", "Ref", "Alt"
    print("\t".join(outheader), file=outhandle)

    handle = open(genotype_count_file, encoding="utf-8")
    x = handle.readline()
    header = x.rstrip("\r\n").split("\t")
    assert header[0] == "Site"
    for x in handle:
        cols = x.rstrip("\r\n").split("\t")
        site = cols[0]
        x = cols[1:]
        counts = [int(x) for x in x]
        if not counts:
            continue
        total = sum(counts)
        if not total:
            continue
        percs = [float(x)/total*100 for x in counts]
        # Nothing can be greater than 100-perc_mixed_calls.
        max_perc = max(percs)
        if max_perc > 100-perc_mixed_calls:
            continue
        # Keep this site.
        # Some chromosome names contain underscores:
        # chr11_KI270831v1_alt_193437_G_A
        x = site.rsplit("_", 3)
        assert len(x) == 4, f"I cannot parse site: {site!r}"
        assert len(x) == len(outheader)
        print("\t".join(x), file=outhandle)


def main():
    genotype_count_file = snakemake.input[0]
    out_file = snakemake.output[0]

    x = int(snakemake.params.keep_sites_with_mixed_calls)
    assert x is None or (x > 0 and x < 100)
    perc_mixed_calls = x

    if perc_mixed_calls is None:
        # Do not filter.  Just list all the sites.
        inhandle = open(genotype_count_file, encoding="utf-8")
        outhandle = open(out_file, 'w', encoding="utf-8")
        for x in inhandle:
            cols = x.rstrip("\r\n").split("\t")
            x = cols[:4]
            print("\t".join(x), file=outhandle)
        inhandle.close()
        outhandle.close()
        return

    # Read in the genotype count file.
    list_good_sites(genotype_count_file, perc_mixed_calls, out_file)


if __name__ == '__main__':
    main()
