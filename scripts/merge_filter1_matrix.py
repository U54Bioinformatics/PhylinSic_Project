def main():
    import os

    matrix_files = snakemake.input
    out_file = snakemake.output[0]

    # Merge the files.  Append the lines, making sure header is
    # the same for every file.
    outhandle = open(out_file, 'w')
    header = None
    for f in matrix_files:
        handle = open(f)
        x = handle.readline()
        assert x
        if not header:
            header = x
            outhandle.write(x)
        assert x == header
        for line in handle:
            outhandle.write(line)
    outhandle.close()


if __name__ == '__main__':
    main()

