def split_list(els, nchunks):

    nel = len(els)

    chunksize = nel // nchunks
    extra_items = nel % nchunks

    chunks = []

    start = 0
    for i in range(nchunks):

        this_chunksize = chunksize
        if i < extra_items:
            this_chunksize += 1

        end = start + this_chunksize

        chunk = els[start:end]
        chunks.append(chunk)

        start = start + this_chunksize

    return [
        chunk for chunk in chunks if len(chunk) > 0
    ]


def write_chunks(chunks, prefix, suffix):
    nchunks = len(chunks)

    cformat = '%0' + str(len(str(nchunks))) + 'i'
    name_format = prefix + cformat + suffix

    for i, chunk in enumerate(chunks):
        fname = name_format % i

        print(fname)

        with open(fname, 'w') as fobj:
            for el in chunk:
                fobj.write(el)
