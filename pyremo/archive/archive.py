


def cdo_call(options='', op='', input='', output=None):
    call = "cdo {} {} {} {}".format(options, op, input, output)
    subprocess.Popen(call, shell=True, stdout=subprocess.PIPE).stdout.read()
    return output
    #return subprocess.run("cdo {} {} {} {}".format(options, op, input, output), shell=True)

    
def convert_with_cdo(f):
    """Convert a single file into NetCDF format.
    """
    file = os.path.basename(f)
    path = os.path.dirname(f)
    #cdo = Cdo()
    #return cdo.copy(options='-f nc', input=f, output=os.path.join(path,'nc',file+'.nc'))
    return cdo_call(options='-f nc', op='copy', input=f, output=os.path.join(path,'nc',file+'.nc'))


def convert_files(files, with_dask=False):
    """Convert many files into NetCDF format.
    """
    results = []
    if with_dask:
        from dask import delayed
    else:
        def delayed(x):
            return x
    for f in files:
        if os.path.isfile(f):
            nc = delayed(convert_with_cdo)(f)
            results.append(nc)
    return results


def get_filename_from_archive(tar, pattern):
    return next(f for f in tar.getnames() if pattern in f)


def extract_file(filename, pattern, path):
    open_tar = tarfile.open(filename, 'r')
    filename = get_filename_from_archive(open_tar, pattern)
    print(filename)
    #dest = os.path.join(scratch, filename)
    open_tar.extract(filename, path=path)
    return os.path.join(path, filename)


def extract_files(tars, pattern="", path=None, parallel=False):
    if path is None:
        path = os.path.join(os.environ['SCRATCH'], '_archive_extract')
    else:
        path = ""
    filenames = []
    if parallel is True:
        from dask import delayed
        futures = []
    else:
        def delayed(x):
            return x
        futures = None
    for tar in tars:
        filenames.append(delayed(extract_file)(tar, pattern, path))
    return filenames