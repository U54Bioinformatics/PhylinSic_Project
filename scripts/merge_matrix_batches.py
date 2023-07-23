
def main():
    import os

    matrix_files = snakemake.input
    out_file = snakemake.output[0]

    # Merge the files.  Append the lines, making sure header is
    # the same for every file.
    outhandle = open(out_file, 'w')
    headers = []  # first 3 lines should be header.
    for f in matrix_files:
        handle = open(f)
        line0 = handle.readline()
        line1 = handle.readline()
        line2 = handle.readline()
        assert line0 and line1 and line2
        if not headers:
            headers = [line0, line1, line2]
            outhandle.writelines(headers)
        assert len(headers) == 3
        assert headers[0] == line0
        assert headers[1] == line1
        assert headers[2] == line2
        for line in handle:
            outhandle.write(line)
    outhandle.close()


if __name__ == '__main__':
    main()

