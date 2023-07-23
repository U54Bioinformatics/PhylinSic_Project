def extract_one(in_file, sample, out_file):
    outhandle = open(out_file, 'w', encoding="utf-8")
    header = "Chrom", "Pos", "Ref", "Alt", "Sample", "Coverage"
    print("\t".join(header), file=outhandle)

    for line in open(in_file, encoding="utf-8"):
        if line.startswith("#"):
            continue

        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) == 10

        chrom, pos = cols[:2]
        ref, alt = cols[3:5]
        assert alt.find(",") < 0, \
               f"Unhandled multiple alt alleles: {line.strip()}"

        # GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
        # ./.:.:1
        assert cols[8] == \
               "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"
        names = cols[8].split(":")
        values = cols[9].split(":")
        if len(values) == 3:   # ignore these sites.  No data.
            continue
        assert len(names) == len(values), f"{cols[8]}\n{cols[9]}"
        RD = values[names.index("RD")]
        AD = values[names.index("AD")]
        FREQ = values[names.index("FREQ")]
        num_ref = int(RD)
        num_alt = int(AD)

        # See a messed up unicode file for FREQ in one of the files.
        # If this happens, then just try to calculate the vaf.
        x = FREQ
        assert x.endswith("%")
        x = x[:-1]
        try:
            vaf = float(x)/100
        except ValueError as x:
            if str(x).startswith("could not convert string to float"):
                assert num_ref + num_alt
                vaf = float(num_alt) / (num_ref + num_alt)
            else:
                raise

        call = f"{num_ref}/{num_alt}/{vaf}"
        x = chrom, pos, ref, alt, sample, call
        assert len(x) == len(header)
        print("\t".join(map(str, x)), file=outhandle)

    outhandle.close()

    
def main():
    import os
    
    in_dir = snakemake.input[0]
    out_dir = snakemake.output[0]
    sample = snakemake.params.sample

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    x = os.listdir(in_dir)
    x1 = [x for x in x if x.endswith(".vcf")]
    x = [os.path.join(in_dir, x) for x in x1]
    in_files = x
    x = [x.replace(".vcf", ".txt") for x in x1]
    x = [os.path.join(out_dir, x) for x in x]
    out_files = x

    for i in range(len(in_files)):
        extract_one(in_files[i], sample, out_files[i])


if __name__ == '__main__':
    main()
