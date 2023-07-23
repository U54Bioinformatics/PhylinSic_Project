def count_genotypes(genotype_file, outfile):
    inhandle = open(genotype_file, encoding="utf-8")
    x = inhandle.readline()
    header = x.rstrip("\r\n").split("\t")
    assert header[0] in ["Site", "Mutation"], \
           f"Unexpected header:{header[0]}\n{genotype_file}"
    site2genotype2count = {}
    for x in inhandle:
        cols = x.rstrip("\r\n").split("\t")
        site = cols[0]
        genotypes = cols[1:]
        genotype2count = {}
        for x in genotypes:
            genotype2count[x] = genotype2count.get(x, 0) + 1
        assert site not in site2genotype2count
        site2genotype2count[site] = genotype2count
    inhandle.close()

    all_genotypes = set()
    for genotype2count in site2genotype2count.values():
        all_genotypes.update(genotype2count)
    all_genotypes = sorted(all_genotypes)

    # Write out the files.
    handle = open(outfile, 'w', encoding="utf-8")
    header = ["Site"] + all_genotypes
    print("\t".join(header), file=handle)
    for site in sorted(site2genotype2count):
        genotype2count = site2genotype2count[site]
        counts = [genotype2count.get(x, 0) for x in all_genotypes]
        x = [site] + counts
        assert len(x) == len(header)
        print("\t".join(map(str, x)), file=handle)
    handle.close()


def main():
    genotype_file = snakemake.input[0]
    genotype_count_file = snakemake.output[0]

    count_genotypes(genotype_file, genotype_count_file)


if __name__ == '__main__':
    main()
