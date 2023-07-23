def main():
    import os

    neighbor_files = snakemake.input
    out_file = snakemake.output[0]
    K = snakemake.params.K

    # Find the nearest neighbors.  For each cell, keep track of the K
    # neighbors with the highest scores.
    MAX_LIST_SIZE = max(K*2, 64)
    scores = {}  # cell1 -> list of (-score, cell)
    for f in neighbor_files:
        it = open(f, encoding="utf-8")

        x = next(it)
        header = x.rstrip("\r\n").split("\t")
        # Header for each file is:
        # Score  [cell name, ...]
        assert header[0] == "Score"
        for x in it:
            cols = x.rstrip("\r\n").split("\t")
            cell1 = cols[0]
            for i in range(1, len(cols)):
                cell2 = header[i]
                if cell1 == cell2:
                    continue
                score = float(cols[i])
                if cell1 not in scores:
                    scores[cell1] = []
                scores[cell1].append((-score, cell2))
                if cell2 not in scores:
                    scores[cell2] = []
                scores[cell2].append((-score, cell1))
            # Keep just the K highest scores.
            for cell in scores:
                if len(scores[cell]) < MAX_LIST_SIZE:
                    continue
                x = sorted(scores[cell])
                scores[cell] = x[:K]

    # Write out the neighbors file.
    handle = open(out_file, 'w', encoding="utf-8")
    x = ["N%02d" % (i+1) for i in range(K)]
    header = ["Cell"] + x
    print("\t".join(header), file=handle)
    for cell in sorted(scores):
        x = sorted(scores[cell])
        x = [x[-1] for x in x]
        x = x[:K]
        assert len(x) == K
        x = [cell] + x
        assert len(x) == len(header)
        print("\t".join(map(str, x)), file=handle)



if __name__ == '__main__':
    main()

