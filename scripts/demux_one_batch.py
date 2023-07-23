# Functions:
# parse_alignment
# get_handle


class SAMAlignment(object):
    # qname    identical for mates.  flag is different.
    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tag_cols):
        # tags should be dict of tag -> (type, value)
        # e.g.:
        # "NH" : ("i", "4")
        assert type(flag) is type(0)
        assert type(pos) is type(0)
        assert type(mapq) is type(0)
        assert type(tlen) is type(0)
        #assert type(tags) is type({})
        assert type(tag_cols) is type([])
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        
        # Optimization: Most of the time, don't care about the tags,
        # but it takes a (relatively) long time to parse.  Don't parse
        # it until needed.
        #self.tags = tags.copy()
        self._tag_cols = tag_cols
        self._tags_parsed = None
        
    def __getattr__(self, attr):
        if attr != "tags":
            raise AttributeError(attr)
        if self._tags_parsed is None:
            tags = {}
            x = [x.split(":", 2) for x in self._tag_cols]
            for (tag, type_, value) in x:
                if tag in tags:
                    # Tophat generates XS tags:
                    # XS:A:-
                    # XS:A:+
                    # Sometimes the same sequence can have both - and
                    # + tags.  If this happens, then just set value to
                    # +-.
                    if tag == "XS" and tags[tag][0] == "A":
                        x = [tags[tag][1], value]
                        x = sorted(set(x))
                        value = "".join(x)
                        tags[tag] = (type_, value)
                        continue
                    #assert tag not in tags, "Duplicate tag %s: %s" % (
                    #    tag, line.rstrip())
                    assert tag not in tags, "Duplicate tag %s: %s" % (
                        tag, self.qname)
                tags[tag] = (type_, value)
            self._tags_parsed = tags
        return self._tags_parsed


def parse_alignment(line, filename=None):
    # Parses a non-header line and return a SAMAlignment.
    # filename is an optional parameter with the name of the bam file.
    # It is only used for error messages.

    # This function is highly optimized, since it may be called many
    # times (e.g. when reading BAM files).
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) < 11:
        # Ignore message:
        # [E::idx_find_and_load] Could not retrieve index file for '/file0/IPCT-S4014-MOON0051-Cap2498-4-NT84_ngs-pipe-2-CACCTTAC.bwa_recalibed.bam'
        
        # Optimization: Don't create the message unless there is a
        # problem.
        msg = "Invalid line (%d):\n%s" % (len(cols), line)
        if filename is not None:
            msg = "Invalid line (%d) from file %s:\n%s" % (
                len(cols), filename, line)
        assert len(cols) >= 11, msg
    qname = cols[0]
    flag = int(cols[1], 10)  # faster conversion if specify base
    rname = cols[2]
    pos = int(cols[3], 10)
    mapq = int(cols[4], 10)
    cigar = cols[5]
    rnext = cols[6]
    pnext = cols[7]
    tlen = int(cols[8], 10)
    seq = cols[9]
    qual = cols[10]
    tag_cols = cols[11:]
    
    # Optimization: create the object myself.  Much faster.
    #x = SAMAlignment(
    #    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq,
    #    qual, tag_cols)
    x = SAMAlignment.__new__(SAMAlignment)
    x.qname = qname
    x.flag = flag
    x.rname = rname
    x.pos = pos
    x.mapq = mapq
    x.cigar = cigar
    x.rnext = rnext
    x.pnext = pnext
    x.tlen = tlen
    x.seq = seq
    x.qual = qual
    x._tag_cols = tag_cols
    x._tags_parsed = None
    return x


def get_handle(barcode, sample, directory, header_lines, barcode2handle):
    # Return a handle to the file.  If the file does not exist yet,
    # create it.  barcode2handle is a dictionary of barcode -> file
    # handle.  Will modify this variable in place.
    import os
    
    if barcode not in barcode2handle:
        x = barcode
        filename = os.path.join(directory, "%s_%s.sam" % (sample, x))
        handle = open(filename, 'w')
        handle.writelines(header_lines)
        barcode2handle[barcode] = handle
    return barcode2handle[barcode]


def main():
    import os
    
    align_file = snakemake.input.align_file
    header_file = snakemake.input.header_file
    out_dir = snakemake.output[0]
    sample = snakemake.params.sample
    
    # If there's 8 batches and 32 handles, can open files for 256
    # cells.  There may be 1000-5000 cells for each 10x run.  Chances
    # of cache miss are pretty high, if cells accessed randomly.
    # In practice, cache misses are pretty rare (< 1:1000).
    MAX_HANDLES = 32

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Read out the SAM file header lines.
    header_lines = open(header_file).readlines()

    barcode2handle = {}
    for line in open(align_file):
        align = parse_alignment(line)
        
        # CR  Cell barcode
        # CY  Cell barcode read quality
        # CB  Cell barcode, error-corrected, confirmed.  
        #     <barcode>-<gem_chip_channel>
        # UR  UMI
        # UY  UMI read quality
        # UB  UMI, error-corrected

        # According to 10x, may be BX:Z for longranger.
        # CB:Z:CCACGAGCATCAACCA-1
        x = align.tags.get("CB", None)
        if not x:   # no barcode
            # Ignore reads without a barcode.
            continue
        x, barcode = x
        assert x == "Z"

        handle = get_handle(
            barcode, sample, out_dir, header_lines, barcode2handle)
        handle.write(line)
        assert len(barcode2handle) <= MAX_HANDLES

    # Close the file handles.
    for handle in barcode2handle.values():
        handle.close()


if __name__ == '__main__':
    main()
