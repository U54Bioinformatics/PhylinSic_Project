class Call:
    # num_ref   int or None
    # num_alt   int or None
    # vaf       float or None
    def __init__(self, num_ref, num_alt, vaf):
        assert num_ref is None or type(num_ref) is type(0)
        assert num_alt is None or type(num_alt) is type(0)
        assert vaf is None or type(vaf) is type(0.0)
        assert num_ref is not None or num_alt is not None or vaf is not None

        total = None
        if num_ref is not None and num_alt is not None:
            total = num_ref + num_alt
        self.num_ref = num_ref
        self.num_alt = num_alt
        self.vaf = vaf
        self.total = total


def _format_call(call):
    # <ref>/<alt>/<vaf>
    # Some callers do not provide VAF.
    if not call:
        return ""
    ref = alt = vaf = ""
    if call.num_ref is not None:
        ref = str(call.num_ref)
    if call.num_alt is not None:
        alt = str(call.num_alt)
    if call.vaf is not None:
        vaf = "%.3f" % call.vaf
    return "%s/%s/%s" % (ref, alt, vaf)


def make_matrix(start, skip, filename, metadata_file, outfile):
    # Read the metadata_file.
    #print "Reading metadata", metadata_file
    coord_data = set()
    i_coord = 0
    samples = set()
    caller2sources = {}
    for x in open(metadata_file):
        cols = x.rstrip("\r\n").split("\t")
        if cols[0] == "Coord":
            chrom, pos, ref, alt = cols[1:]
            pos = int(pos)
            if i_coord % skip == start:
                coord_data.add((chrom, pos, ref, alt))
            i_coord += 1
        elif cols[0] == "Source":
            caller, source = cols[1:]
            if caller not in caller2sources:
                caller2sources[caller] = set()
            caller2sources[caller].add(source)
        elif cols[0] == "Sample":
            samples.add(cols[1])
        else:
            raise AssertionError("I don't understand row: %s" % repr(cols))
    samples = sorted(samples)
    callers = sorted(caller2sources)
    #print("Found %d coordinates, %d samples, %d callers" % (
    #    len(coord_data), len(samples), len(callers)))

    # Make a list of all calls.
    #print "Reading calls", filename
    call_data = []
    it = open(filename)
    header = next(it).rstrip("\r\n").split("\t")
    i_Chrom = header.index("Chrom")
    i_Pos = header.index("Pos")
    i_Ref = header.index("Ref")
    i_Alt = header.index("Alt")
    i_Num_Alt = header.index("Num Alt")
    i_Num_Ref = header.index("Num Ref")
    i_VAF = header.index("VAF")
    i_Sample = header.index("Sample")
    i_Caller = header.index("Caller")
    i_Source = header.index("Source")
    for x in it:
        cols = x.rstrip("\r\n").split("\t")
        x = (cols[i_Chrom], int(cols[i_Pos]), cols[i_Ref], cols[i_Alt])
        if x not in coord_data:
            continue

        assert cols[i_Num_Alt].find(",") < 0

        # Get the calls.
        assert cols[i_Source] in ["DNA", "RNA"], "Unknown source: %s" % \
               cols[i_Source]
        num_ref = num_alt = vaf = None
        if cols[i_Num_Ref]:
            num_ref = int(cols[i_Num_Ref])
        if cols[i_Num_Alt]:
            num_alt = int(cols[i_Num_Alt])
        if cols[i_VAF]:
            vaf = float(cols[i_VAF])
        if num_ref is None and num_alt is None and vaf is None:
            continue
        call = Call(num_ref, num_alt, vaf)
        x = (cols[i_Chrom], int(cols[i_Pos]), cols[i_Ref], cols[i_Alt],
             cols[i_Sample], cols[i_Caller], cols[i_Source], call)
        call_data.append(x)

    #print "Processing %d calls" % len(call_data)
    # sample -> caller -> source -> chrom, pos, ref, alt -> call
    samp2caller2source2coord2call = {}
    for x in call_data:
        chrom, pos, ref, alt, sample, caller, source, call = x
        coord = chrom, pos, ref, alt
        if sample not in samp2caller2source2coord2call:
            samp2caller2source2coord2call[sample] = {}
        caller2source2coord2call = samp2caller2source2coord2call[sample]
        if caller not in caller2source2coord2call:
            caller2source2coord2call[caller] = {}
        source2coord2call = caller2source2coord2call[caller]
        if source not in source2coord2call:
            source2coord2call[source] = {}
        coord2call = source2coord2call[source]

        # A (sample, caller, coord) may have multiple calls.  For
        # example, for germline samples that are called with each
        # tumor sample.  If this is the case, then take the call
        # with the highest coverage.
        if coord in coord2call:
            old_call = coord2call[coord]
            cov = old_cov = None
            if call.num_ref is not None and call.num_alt is not None:
                cov = call.num_ref + call.num_alt
            if old_call.num_ref is not None and \
                   old_call.num_alt is not None:
                old_cov = old_call.num_ref + old_call.num_alt
            if cov is None and old_cov is not None:
                call = old_call
            elif cov is not None and old_cov is not None and cov < old_cov:
                call = old_call
        coord2call[coord] = call

    # Count the number of callers that called a variant at each
    # position for each sample.
    # sample -> chrom, pos, ref, alt -> caller -> source -> 1
    samp2coord2caller2source = {}
    # Need to do this first, to make sure each caller is counted
    # at most once.  This is to account for germline samples that
    # is called by each caller multiple times.
    for x in call_data:
        chrom, pos, ref, alt, sample, caller, source, call = x
        coord = chrom, pos, ref, alt
        if sample not in samp2coord2caller2source:
            samp2coord2caller2source[sample] = {}
        if coord not in samp2coord2caller2source[sample]:
            samp2coord2caller2source[sample][coord] = {}
        if caller not in samp2coord2caller2source[sample][coord]:
            samp2coord2caller2source[sample][coord][caller] = {}
        samp2coord2caller2source[sample][coord][caller][source] = 1
    samp2coord2nc = {}  # sample -> chrom, pos, ref, alt -> num_callers
    for sample in samp2coord2caller2source:
        samp2coord2nc[sample] = {}
        for coord in samp2coord2caller2source[sample]:
            caller2source = samp2coord2caller2source[sample][coord]
            count = 0
            for caller in caller2source:
                count += len(caller2source[caller])
            samp2coord2nc[sample][coord] = count

    #print "Formatting matrix"
    # Format everything into an annotation matrix.
    coord_data = sorted(coord_data)
    headers0 = []
    headers1 = []
    headers2 = []
    all_annots = []

    # Add the positions.
    headers0 += ["", "", "", ""]
    headers1 += ["", "", "", ""]
    headers2 += ["Chrom", "Pos", "Ref", "Alt"]
    for i in range(4):
        x = [x[i] for x in coord_data]
        x = [str(x) for x in x]
        all_annots.append(x)

    # Add the number of callers information.
    headers0 += ["Num Callers"] * len(samples)
    headers1 += [""] * len(samples)
    headers2 += samples
    for sample in samples:
        annots = []
        for coord in coord_data:
            nc = samp2coord2nc.get(sample, {}).get(coord, "")
            annots.append(nc)
        all_annots.append(annots)

    # Add information about calls.
    for sample in samples:
        #print("Getting calls for sample", sample)
        sample_printed = False
        x = samp2caller2source2coord2call.get(sample, {})
        caller2source2coord2call = x
        for caller in callers:
            for source in sorted(caller2sources[caller]):
                h0 = ""
                if not sample_printed:
                    h0 = sample
                    sample_printed = True
                h1 = caller
                if source != "DNA":
                    h1 = "%s (%s)" % (caller, source)
                h2 = "Ref/Alt/VAF"
                headers0.append(h0)
                headers1.append(h1)
                headers2.append(h2)

                x = caller2source2coord2call.get(caller, {}).get(source, {})
                coord2call = x
                annots = []
                for coord in coord_data:
                    x = ""
                    call = coord2call.get(coord)
                    if call:
                        x = _format_call(call)
                    annots.append(x)

                all_annots.append(annots)


    #print "Writing matrix", outfile
    # Set the headers.
    assert len(headers0) == len(headers1)
    assert len(headers0) == len(headers2)
    assert len(headers0) == len(all_annots)

    for i in range(len(headers0)-1, -1, -1):
        # If headers1[i] is the same as header1[i-1], then do not
        # write it out again.
        #
        # Exception: If headers0[i] != headers0[i-1], then we're
        # starting a new "block", and headers1[i] should still be
        # written out.
        # Example: If there's only one <Caller>, then the <Sample>
        # will not be blank, but the <Caller> should still be copied
        # over (because they are the same).
        # <Sample1>    <Sample2>
        # <Caller>     <Caller>
        # Ref/Alt/VAF  Ref/Alt/VAF
        if headers1[i] == headers1[i-1] and headers0[i] == headers0[i-1]:
            headers1[i] = ""
        if headers0[i] == headers0[i-1]:
            headers0[i] = ""

    handle = open(outfile, 'w')
    print("\t".join(headers0), file=handle)
    print("\t".join(headers1), file=handle)
    print("\t".join(headers2), file=handle)
    for i in range(len(all_annots[0])):
        x = [x[i] for x in all_annots]
        print("\t".join(map(str, x)), file=handle)

    #print "Done"
    


def main():
    variant_file = snakemake.input[0]
    metadata_file = snakemake.input[1]
    out_file = snakemake.output[0]
    batch = int(snakemake.params.batch)
    NUM_BATCHES = snakemake.params.num_batches

    make_matrix(batch, NUM_BATCHES, variant_file, metadata_file, out_file)


if __name__ == '__main__':
    main()

