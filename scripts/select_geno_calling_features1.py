def choose_features(site_meta_file, num_features, outfile):
    # Select the features used to calculate the neighbors.
    MIN_MINOR_CELLS = 5

    inhandle = open(site_meta_file, encoding="utf-8")
    x = inhandle.readline()
    header = x.rstrip("\r\n").split("\t")
    mut_data = []
    for x in inhandle:
        x = x.rstrip("\r\n").split("\t")
        mut_data.append(x)

    # Only keep sites where at least MIN_MINOR_CELLS have the minor
    # genotype.
    # On test data set, MIN_MINOR_CELLS=5 got rid of 94% of sites.
    # MIN_MINOR_CELLS=1 got rid of 73% of sites.
    I_Num_Minor = header.index("Num Minor")
    x = [x for x in mut_data if int(x[I_Num_Minor]) >= MIN_MINOR_CELLS]
    if len(x) >= num_features:
        # Only use this cutoff if there are enough features.
        mut_data = x

    # Choose the sites with the highest number of reads.
    I_Total_Reads = header.index("Total Reads")
    schwartz = [(-int(x[I_Total_Reads]), x) for x in mut_data]
    schwartz.sort()
    mut_data = [x[-1] for x in schwartz]
    mut_data = mut_data[:num_features]

    # Write out the feature file.
    I_Site = header.index("Site")
    x = [x[I_Site] for x in mut_data]
    x = [f"{x}\n" for x in x]
    with open(outfile, 'w', encoding="utf-8") as h:
        h.writelines(x)


def main():
    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    ternary_features = snakemake.params.ternary_features
    assert ternary_features >= 1 and ternary_features <= 1024

    choose_features(in_file, ternary_features, out_file)


if __name__ == '__main__':
    main()
