PARSE_RAV_CACHE = {
    #"1/0/0.000" : (1, 0, 0.0),
    #"2/0/0.000" : (2, 0, 0.0),
    #"3/0/0.000" : (3, 0, 0.0),
    #"0/1/0.000" : (0, 1, 0.0),
    #"0/2/0.000" : (0, 2, 0.0),
    #"0/3/0.000" : (0, 3, 0.0),
    }
def parse_rav(rav_value):
    # rav_value should be <ref>/<alt>/<vaf>.  Return a tuple of <ref>,
    # <alt>, <vaf>.  <ref> and <alt> are ints, and <vaf> is a float.

    # Caching makes it ~12% faster.  Not a huge speedup, but helps for
    # large files.

    # 15.598  No memoization.
    # 14.353  Only save known values.
    # 13.754  Memoization.
    x = PARSE_RAV_CACHE.get(rav_value)
    if x is None:  # function never returns None
        x = _parse_rav_h(rav_value)
        PARSE_RAV_CACHE[rav_value] = x
    return x


def _parse_rav_h(rav_value):
    if not rav_value:
        return 0, 0, 0.0
    x = rav_value.split("/")
    assert len(x) == 3, "Does not look like <ref>/<alt>/<vaf>: %s" % rav_value
    ref, alt, vaf = x
    # Handle case where RAV is: "//0.000"
    if not ref and not alt:
        ref = alt = 0
    try:
        ref, alt, vaf = int(ref), int(alt), float(vaf)
    except ValueError as x:
        raise ValueError(
            "I could not parse <ref>/<alt>/<vaf>: %s\n%s" % (
                rav_value, str(x)))
    return ref, alt, vaf


def main():
    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    # Keep a site if at least one of the cells has this number of
    # reads.
    keep_sites_reads = snakemake.params.keep_sites_with_high_reads
    if keep_sites_reads is not None:
        assert keep_sites_reads >= 1 and keep_sites_reads < 16384

    x = snakemake.params.remove_chrom_with_prefix or ""
    x = x.split(",")
    x = [x for x in x if x]
    discard_prefixes = set(x)

    x = snakemake.params.remove_multiallelic_sites
    assert x in ["yes", "no", None], \
           'remove_multiallelic_sites must be "yes" or "no"'
    remove_multiallelic = x == "yes"

    inhandle = open(in_file, encoding="utf-8")
    outhandle = open(out_file, 'w', encoding="utf-8")

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

        chrom, pos, ref, alt = cols[:4]
        pos = int(pos)

        # Ignore this site if the prefix starts with any of the
        # prefixes in discard_prefixes.
        x = [x for x in discard_prefixes if chrom.startswith(x)]
        if x:
            continue

        if remove_multiallelic:
            if ref.find(",") >= 0:
                continue
            if alt.find(",") >= 0:
                continue

        if keep_sites_reads is not None:
            found = False
            for x in cols[4:]:
                if not x:
                    continue
                r, a, v = parse_rav(x)
                total = r + a

                if total >= keep_sites_reads:
                    found = True
                    break

            if not found:
                continue

        print(line, file=outhandle, end="")

    inhandle.close()
    outhandle.close()




if __name__ == '__main__':
    main()
