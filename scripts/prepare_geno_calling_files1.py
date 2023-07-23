CALL_GENOTYPE_DICT = {}
def call_genotype_ternary(
    ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff, prob_cutoff):
    # Return None, "REF", "ALT", or "HET".  vaf_cutoff determines the
    # variant allele frequency needed to call something heterozygous.
    # For instance, vaf_cutoff of 0.3 means 0-0.3 is REF, 0.3-0.7 is
    # HET, and 0.7-1.0 is ALT.  If the probability of the genotype <
    # prob_cutoff, will return None.  prob_cutoff should be from [0,
    # 100].
    key = ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff, prob_cutoff
    if key not in CALL_GENOTYPE_DICT:
        CALL_GENOTYPE_DICT[key] = _call_genotype_ternary_h(
            ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff,
            prob_cutoff)
    return CALL_GENOTYPE_DICT[key]


def _call_genotype_ternary_h(
    ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff, prob_cutoff):
    assert prob_cutoff > 1 and prob_cutoff <= 100
    scores = score_genotype_ternary(
        ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff)
    best_call = best_p = None
    for call, p in scores.items():
        if best_p is None or p > best_p:
            best_call = call
            best_p = p
    if best_p < prob_cutoff:
        best_call = None
    return best_call


def score_genotype_ternary(
    ref_count, alt_count, ref_prior, alt_prior, vaf_cutoff):
    # Return a dictionary of "REF", "ALT", or "HET" -> probability.
    # The probability goes from 0-100.
    import scipy.stats

    assert vaf_cutoff >= 0 and vaf_cutoff < 0.5

    alpha = ref_count+1+ref_prior
    beta = alt_count+1+alt_prior
    x1 = scipy.stats.beta.cdf(vaf_cutoff, alpha, beta) * 100
    x2 = scipy.stats.beta.cdf(1.0-vaf_cutoff, alpha, beta) * 100
    p_ref = 100 - x2
    p_alt = x1
    p_het = 100 - p_ref - p_alt
    ## Can generate warning (sometimes):
    ## 3806 92 0.25 0.25 0.1
    #p_ref = integrate_beta_pdf(
    #    ref_count+ref_prior, alt_count+alt_prior, 1.0-vaf_cutoff, 1.0) * 100
    #p_alt = integrate_beta_pdf(
    #    ref_count+ref_prior, alt_count+alt_prior, 0, vaf_cutoff) * 100
    #p_het = 100 - p_ref - p_alt
    p_ref = min(max(p_ref, 0), 100)
    p_alt = min(max(p_alt, 0), 100)
    p_het = min(max(p_het, 0), 100)
    x = { "REF" : p_ref, "ALT" : p_alt, "HET" : p_het }
    return x


def prepare_files(svm_file, site_meta_file, ref_matrix_file, alt_matrix_file):
    # site_meta_file
    # - Site
    # - Total Reads    Total number of reads across all cells for this site.
    # - Num Cells      Number of cells with >= 1 read.
    # - Variance       Variance of VAF for num cells.
    # - Calls          Number of cells with confident ternary calls.
    # - Perc Major     Percent of cells with major genotype.
    # - Num Minor      Number of cells with minor genotype.
    import numpy
    import filter_sites1

    # Keep track of cells_per_site.
    site2ncells = {}     # site -> num cells
    # Keep track of total reads.
    site2nreads = {}     # site -> number of reads
    # Keep track of variance of VAF.
    site2varvaf = {}     # site -> variance of VAFs
    # Keep track of homogeneity.
    site2ncalled = {}    # site -> num cells used for homogeneous
    site2percmajor = {}  # site -> percent of most common genotype
    site2nminor = {}     # site -> cells with minor genotypes

    it = open(svm_file, encoding="utf-8")
    x = it.readline()
    header = x.rstrip("\r\n").split("\t")
    assert len(header) > 4
    samples = header[4:]

    ref_handle = open(ref_matrix_file, 'w', encoding="utf-8")
    alt_handle = open(alt_matrix_file, 'w', encoding="utf-8")
    header = ["Site"] + samples
    print("\t".join(header), file=ref_handle)
    print("\t".join(header), file=alt_handle)

    for line in it:
        cols = line.rstrip("\r\n").split("\t")
        x = cols[:4]
        site = "_".join(x)

        ref_counts = [""] * (len(cols)-4)
        alt_counts = [""] * (len(cols)-4)
        genotypes = [""] * (len(cols)-4)
        cov_cols = cols[4:]

        vafs = []
        for i, cov in enumerate(cov_cols):
            if not cov:
                continue
            r, a, v = filter_sites1.parse_rav(cov)
            ref_counts[i] = r
            alt_counts[i] = a

            vafs.append(v)
            # May be None.
            # 80% probability is good cutoff.  Assigns 4/4 as
            # heterozygous.  90% is too strict.
            genotypes[i] = call_genotype_ternary(r, a, 0, 0, 0.3, 80)

        # Write out reference and alt reads.
        x1 = [site] + ref_counts
        x2 = [site] + alt_counts
        assert len(x1) == len(header)
        assert len(x2) == len(header)
        print("\t".join(map(str, x1)), file=ref_handle)
        print("\t".join(map(str, x2)), file=alt_handle)

        # Number of expressed cells.
        x = [x for x in ref_counts if x != ""]
        site2ncells[site] = len(x)

        # Total reads.
        x1 = [x for x in ref_counts if x]
        x2 = [x for x in alt_counts if x]
        site2nreads[site] = sum(x1) + sum(x2)

        # Variance.  Variance of the variant allele frequencies.
        varvaf = 0
        if len(vafs) > 1:
            varvaf = numpy.var(vafs)
        site2varvaf[site] = varvaf

        # Homogeneity.  Percent of the most common genotype.
        genotypes = [x for x in genotypes if x]
        perc_major = 100
        nminor = 0
        if genotypes:
            geno2count = {}
            for x in genotypes:
                geno2count[x] = geno2count.get(x, 0) + 1
            x = [(count, geno) for (geno, count) in geno2count.items()]
            x.sort()
            counts = [x[0] for x in x]
            max_count = counts[-1]
            perc_major = int(float(max_count) / len(genotypes) * 100)
            nminor = sum(counts[:-1])
        site2ncalled[site] = len(genotypes)
        site2percmajor[site] = perc_major
        site2nminor[site] = nminor

    # Write out the metadata file.
    handle = open(site_meta_file, 'w', encoding="utf-8")
    header = "Site", "Total Reads", "Num Cells", "Variance", "Called", \
             "Perc Major", "Num Minor"
    print("\t".join(header), file=handle)
    for site, num_cells in site2ncells.items():
        x = site, site2nreads[site], num_cells, site2varvaf[site], \
            site2ncalled[site], site2percmajor[site], site2nminor[site]
        assert len(x) == len(header)
        print("\t".join(map(str, x)), file=handle)


def main():
    import os
    
    in_file = snakemake.input[0]
    assert len(snakemake.output) == 3
    metadata_file, ref_file, alt_file = snakemake.output

    for out_file in snakemake.output:
        out_dir = os.path.dirname(out_file)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    prepare_files(in_file, metadata_file, ref_file, alt_file)


if __name__ == '__main__':
    main()
