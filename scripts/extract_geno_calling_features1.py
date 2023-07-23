def extract_features_from_matrix(
    matrix_file, feature_file, feature_matrix_file):
    
    # Read in the features.
    x = open(feature_file, encoding="utf-8")
    x = [x.strip() for x in x]
    features = set(x)

    seen = set()

    # Write out the ref matrix.
    inhandle = open(matrix_file, encoding="utf-8")
    outhandle = open(feature_matrix_file, 'w', encoding="utf-8")
    header = next(inhandle).rstrip("\r\n").split("\t")
    print("\t".join(header), file=outhandle)
    for x in inhandle:
        cols = x.rstrip("\r\n").split("\t")
        if cols[0] not in features:
            continue
        seen.add(cols[0])
        print("\t".join(cols), file=outhandle)

    assert features == seen, "Missing"


def main():
    assert len(snakemake.input) == 3
    assert len(snakemake.output) == 2

    ref_matrix_file = snakemake.input[0]
    alt_matrix_file = snakemake.input[1]
    feature_file = snakemake.input[2]
    ref_feature_matrix_file = snakemake.output[0]
    alt_feature_matrix_file = snakemake.output[1]
    
    extract_features_from_matrix(
        ref_matrix_file, feature_file, ref_feature_matrix_file)
    extract_features_from_matrix(
        alt_matrix_file, feature_file, alt_feature_matrix_file)


if __name__ == '__main__':
    main()
