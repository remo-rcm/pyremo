


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