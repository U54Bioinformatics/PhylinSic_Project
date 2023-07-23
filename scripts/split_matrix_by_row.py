def main():
    import os
    
    in_file = snakemake.input[0]
    out_files = snakemake.output[:]

    # For large files, we'll run into the filehandle limit.  Do only
    # one batch at a time.
    max_open_files = 32

    for out_file in out_files:
        out_dir = os.path.dirname(out_file)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    num_files = len(out_files)
    num_batches = int(num_files/max_open_files) + 1
    for i in range(num_batches):
        out_handles = [None] * len(out_files)
        for j in range(len(out_files)):
            if j % num_batches == i:
                out_handles[j] = open(out_files[j], 'w', encoding="utf-8")

        # Open the file.
        handle = open(in_file, encoding="utf-8")
        # Write the header lines.
        x = handle.readline()
        assert x
        for h in out_handles:
            if h is not None:
                h.write(x)
        # Write the contents.
        for i, x in enumerate(handle):
            j = i % len(out_handles)
            if out_handles[j] is not None:
                out_handles[j].write(x)

        handle.close()
        for x in out_handles:
            if x is not None:
                x.close()


if __name__ == '__main__':
    main()
