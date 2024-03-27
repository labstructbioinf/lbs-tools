import tempfile, gzip

def get_unzipped_tempfile(gz_fn):
    f = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
    fh = open(f.name, 'w')
    fh.write(gzip.open(gz_fn, mode='rt').read())
    fh.close()
    return f.name

