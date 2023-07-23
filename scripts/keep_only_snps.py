def main():
    import os
    
    in_file = snakemake.input[0]
    out_file = snakemake.output[0]

    handle = open(out_file, 'w', encoding="utf-8")
    for line in open(in_file):
        if not line.strip():
            continue
        if line.startswith("#"):
            handle.write(line)
            continue
        cols = line.rstrip("\r\n").split("\t")
        # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  015_R_LIVER
        assert len(cols) >= 9
        ref = cols[3]
        alt = cols[4]
        # No Indels
        #
        # REF   ALT
        #  C     T     snp
        #  TC    T     indel
        if len(ref) > 1 or len(alt) > 1:
            continue
        handle.write(line)
    handle.close()
        
    
if __name__ == '__main__':
    main()
