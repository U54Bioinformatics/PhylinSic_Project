def split_file(in_file, file_num, num_little_files, out_file):
    # Split a big SVM file into a bunch of little ones, each with a
    # subset of the variants.
    # num_little_files is the number of little files this should be
    # split into.  file_num is the which little file we are creating
    # with this function call.  Should go from 0-(num_little_files-1)
    assert file_num >= 0 and file_num < num_little_files

    inhandle = open(in_file, encoding="utf-8")
    outhandle = open(out_file, 'w', encoding="utf-8")

    line0 = inhandle.readline()
    assert line0
    outhandle.write(line0)

    for i, line in enumerate(inhandle):
        if i % num_little_files != file_num:
            continue
        outhandle.write(line)
    inhandle.close()  # so gzip doesn't live forever
    outhandle.close()


def main():
    import subprocess

    in_file = snakemake.input[0]
    out_file = snakemake.output[0]
    batch = int(snakemake.params.batch)
    NUM_BATCHES = snakemake.params.num_batches

    split_file(in_file, batch, NUM_BATCHES, out_file)


if __name__ == '__main__':
    main()
