=====
Usage
=====

The pyremo package contains some common tools for working with `REMO <www.remo-rcm.de>`_ data.

This includes:

* accessing meta data and tables online,
* working with REMO domains,
* basic REMO specific physics and pressure interpolation,
* and a new preprocessing workflow.

Most of the packages API relies on `xarray data structures
<http://xarray.pydata.org/en/stable/user-guide/data-structures.html>`_.
This means that most of the functions require an xarray data structure
as input and will return xarray data structures. This is very convenient
since xarray will automatically handle dimensions correctly. Xarray will
also automatically handle vectorizations (e.g. do all computations automatically
along the time axis) and also `works very well with dask
<http://xarray.pydata.org/en/stable/user-guide/dask.html>`_ for parallel and out of core
computations.

Most examples and tutorials are available as `notebooks
<https://nbviewer.jupyter.org/github/remo-rcm/pyremo/tree/master/notebooks/>`_
and are also rendered as part of the documentation for convenience.
