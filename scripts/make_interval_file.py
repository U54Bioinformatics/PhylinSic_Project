def read_fa_index(filename):
    # Iterate over each line in the index.
    import os

    for line in open(filename):
        x = line.rstrip("\r\n").split("\t")
        assert len(x) == 5
        name, length, offset, linebases, linewidth = x
        x = type("", (), {})()
        x.name = name
        x.length = int(length)
        x.offset = int(offset)
        x.linebases = int(linebases)
        x.linewidth = int(linewidth)
        yield x


def main():
    import os
    import math
    
    genome_index_file = snakemake.input[0]
    out_file = snakemake.output[0]
    interval = int(snakemake.params.interval)
    num_intervals = snakemake.params.num_intervals
    assert interval >= 0 and interval < num_intervals

    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Coordinates here are 1-based, inclusive end.

    # Read the index.
    data = []
    for x in read_fa_index(genome_index_file):
        d = type("", (), {})()
        d.chrom = x.name
        d.start = 1
        d.end = x.length
        data.append(d)

    # Calculate the number of bases per interval.
    genome_size = sum(x.end for x in data)
    bases_per_interval = math.ceil(genome_size/num_intervals)

    # Split into multiple files of bases_per_interval bases.
    filenum = 0
    filenum2intervals = {}
    intervals = []  # list of (chrom, start, end)
    while data:

        # Calculate the base pairs used in this set of intervals.  If
        # long enough, then split this as a separate set of intervals.
        x = [x[2]-x[1]+1 for x in intervals]
        interval_len = sum(x)
        if interval_len >= bases_per_interval:
            filenum2intervals[filenum] = intervals
            filenum += 1
            intervals = []
            continue

        # Look at the next chromosome.
        d = data.pop(0)

        # Calculate the length of this chromosome (or part of chromosome).
        chrom_len = d.end - d.start + 1

        if interval_len + chrom_len <= bases_per_interval:
            # If this whole (rest of this) chromosome can fit within
            # this interval, then add the rest of this chromosome.
            intervals.append((d.chrom, d.start, d.end))
        else:
            # Add a piece of this chromosome.
            bases_to_add = bases_per_interval - interval_len
            assert bases_to_add > 0
            end = d.start + bases_to_add - 1
            assert end < d.end
            intervals.append((d.chrom, d.start, end))

            # Set the new start to the previous end.
            x = type("", (), {})()
            x.chrom = d.chrom
            x.start = end + 1
            x.end = d.end
            data.insert(0, x)
    if intervals:
        filenum2intervals[filenum] = intervals

    assert len(filenum2intervals) == num_intervals, "%d %d" % (
        len(filenum2intervals), num_intervals)
    intervals = filenum2intervals[interval]

    handle = open(out_file, 'w')
    for (chrom, start, end) in intervals:
        print("%s:%d-%d" % (chrom, start, end), file=handle)
    handle.close()
        
    
if __name__ == '__main__':
    main()
