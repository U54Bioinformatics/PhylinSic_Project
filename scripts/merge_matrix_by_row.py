def safe_unlink(filename):
    # Will remove a file or directory.
    import os
    import shutil

    if not filename:
        return
    if not os.path.islink(filename) and not os.path.exists(filename):
        return
    try:
        if os.path.isfile(filename) or os.path.islink(filename):
            os.unlink(filename)
        else:
            shutil.rmtree(filename)
    except (OSError, IOError) as x:
        if str(x).find("No such file") >= 0:
            pass
        else:
            raise


class Temphandle:
    def __init__(self, suffix=None, prefix=None, dir=None, unlink=True):
        import os
        import tempfile
        params = {}
        if suffix is not None:
            params["suffix"] = suffix
        if prefix is not None:
            params["prefix"] = prefix
        if dir is not None:
            params["dir"] = dir
        x, filename = tempfile.mkstemp(**params)
        os.close(x)
        self.filename = filename
        self.unlink = unlink
    def __enter__(self):
        return self.filename
    def __exit__(self, exc_type, exc_value, tb):
        import traceback
        if self.unlink:
            safe_unlink(self.filename)
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)
            return False
        return True
    def __del__(self):
        import os
        # TODO: Fix error:
        if self.unlink:
            # No.  This generates an error for some reason.
            # Exception RuntimeError: 'sys.path must be a list of
            # directory names' in <bound method Tempfile.__del__ of
            # <genomicode.Tempfile.Tempfile instance at 0x2aaab8510e18>>
            # ignored
            try:
                os.unlink(self.filename)
            except (OSError, IOError) as x:
                if str(x).find("No such file") >= 0:
                    pass
                else:
                    raise


# split_by_row and merge_by_row not guaranteed to maintain the order
# of the rows.
def merge_by_row(filenames, outfile, num_header_rows=None,
                 max_open_files=None):
    import shutil

    if not filenames:
        return
    if len(filenames) == 1:
        shutil.copy2(filenames[0], outfile)
        return

    hr = num_header_rows
    assert hr

    if max_open_files is None:
        max_open_files = 64
    assert max_open_files is None or max_open_files >= 1

    if len(filenames) <= max_open_files:
        _merge_by_row_h(filenames, outfile, hr)
        return

    # If there are more files than the file handle limit, then merge
    # these files recursively.
    # Split filenames into batches of max_open_files.
    num_batches = len(filenames) / max_open_files + 1
    assert num_batches >= 1

    batch_filenames = []  # list of lists
    batch_outhandles = []
    for i in range(num_batches):
        start = i * max_open_files
        end = start + max_open_files
        x = filenames[start:end]
        batch_filenames.append(x)
        x = Temphandle(suffix=".txt", dir=".")
        batch_outhandles.append(x)
    # Merge each of the batches.
    for i in range(num_batches):
        # Can use _merge_by_row_h here because we know each batch has
        # less than max_open_files files.
        _merge_by_row_h(batch_filenames[i], batch_outhandles[i].filename, hr)

    # Merge the batch outfiles to the final outfile.  Since there may
    # be more batches than max_open_files, call merge_by_row
    # recursively.
    x = [x.filename for x in batch_outhandles]
    merge_by_row(
        x, outfile, num_header_rows=hr, max_open_files=max_open_files)


def _merge_by_row_h(filenames, outfile, num_header_rows):
    outhandle = open(outfile, 'w', encoding="utf-8")
    handles = [open(x, encoding="utf-8") for x in filenames]

    # Read the header lines and make sure they are the same in every
    # file.
    for i in range(num_header_rows):
        lines = [x.readline() for x in handles]
        for j in range(1, len(lines)):
            assert lines[0] == lines[j]
        outhandle.write(lines[0])

    while handles:
        i = 0
        while i < len(handles):
            x = handles[i].readline()
            if not x:
                del handles[i]
                continue
            outhandle.write(x)
            i += 1

    for x in handles:
        x.close()
    outhandle.close()



def main():
    in_files = snakemake.input[:]
    out_file = snakemake.output[0]

    merge_by_row(in_files, out_file, num_header_rows=1, max_open_files=32)


if __name__ == '__main__':
    main()
