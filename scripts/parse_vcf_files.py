def main():
    import os

    in_files = snakemake.input
    out_file = snakemake.output[0]

    handle = open(out_file, 'w', encoding="utf-8")

    header = (
        "Caller", "File", "Sample", "Chrom", "Pos", "Ref", "Alt", "Source",
        "Num Ref", "Num Alt", "Total Reads", "VAF", "Filter", "Call", "GQ")
    print("\t".join(header), file=handle)

    caller_name = "HaplotypeCaller"
    source = "RNA"

    for vcf_file in in_files:
        filestem = os.path.split(vcf_file)[1]
        sample = None
        for line in open(vcf_file, encoding="utf-8"):
            if line.startswith("#CHROM"):
                cols = line.rstrip("\r\n").split("\t")
                assert len(cols) == 10
                sample = cols[9]
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  015_R_LIVER
            assert len(cols) == 10
            chrom = cols[0]
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            filter_str = cols[6]

            # INFO
            # AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=20.83;SOR=4.174
            # FORMAT
            # GT:AD:DP:GQ:PL
            # 015_R_LIVER
            # 1/1:0,7:7:21:174,21,0

            format_spec = cols[8]
            assert format_spec == "GT:AD:DP:GQ:PL"
            x = cols[9]
            x = x.split(":")
            assert len(x) == 5
            call_str, AD, total_reads, GQ, x = x
            x = AD.split(",")
            assert len(x) == 2
            num_ref, num_alt = int(x[0]), int(x[1])
            vaf = 0
            if num_alt + num_ref:
                vaf = float(num_alt) / (num_alt + num_ref)

            assert sample is not None
            x = caller_name, filestem, sample, chrom, pos, \
                ref, alt, source, \
                num_ref, num_alt, total_reads, vaf, filter_str, call_str, GQ
            assert len(x) == len(header)
            x = "\t".join(map(str, x))
            print(x, file=handle)



if __name__ == '__main__':
    main()
