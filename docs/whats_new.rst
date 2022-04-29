.. currentmodule:: pyremo

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo

v0.3.0 (Unreleased)
-------------------------

New Features
~~~~~~~~~~~~

- ``prsint``: a command line interface for the pressure interpolator (:pull:`53`).
- ERA5 cmorizer includes ERA5.1 to work with DKRZ data pool (:pull:`52`).
- Added ``python3.10`` support (:pull:`46`).
- Updated cmorization module for use with CMIP6 (:pull:`48`, :pull:`49`, :pull:`51`) and CORDEX vocabulary. The underlying tables are only used for testing and should not yet be uses for actual data publication, see also `here <https://github.com/euro-cordex/py-cordex/pull/55>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Updated documentation with improved contribution guide (:pull:`63`).
- Fixed issues with ``dask.delayed`` in ``preproc.era5`` running in batch mode.
- Added ``.pre-commit-config.yaml`` and ``linting.yaml`` to run linter checks (:pull:`55`).
- Code is reformatted to apply to ``flake8`` and ``black`` conventions (:pull:`55`).
- ERA5 cmorizer works with pandas datatables (:pull:`52`).
- Tutorial `data source <https://github.com/remo-rcm/pyremo-data>`_ is now on github (:pull:`47`).
- Added fortran extensions ``pyintorg`` and ``pydruint`` to extra dependencies.
- Added ``CI-extensions.yaml`` for git actions testing with fortran extensions.
- Included tests for ``pyremo.prsint`` with fortran extension (:pull:`43`).
- Pinned ``sphinx`` and ``jinja`` dependencies (:pull:`45`).

Breaking Changes
~~~~~~~~~~~~~~~~

- Updated cmor API of ``cmor.cmorize_variable`` to use actual filenames of tables (:pull:`48`).

Bugfixes
~~~~~~~~

- Fixed bug in ``prsint`` for vertical level coordinate (:pull:`57`).


v0.2.0 (24 February 2022)
-------------------------

This is a major restructuring release and comes with a new package structure and a lot of reintegrated tools.
It includes a new preprocessing module, pressure interpolation and cmorization based on an xarray API.

New Features
~~~~~~~~~~~~
- Included physics package with xarray API.
- Included preprocessing interface for xarray data structures.
- Included ERA5 cmorizer.
- Includes experimental cmorization module (:pull:`20`, :pull:`21`, :pull:`22`,
  :pull:`23`, :pull:`25`, :pull:`26`).
- Included production analysis (:pull:`33`).
- Command line tool for remo analysis (:pull:`35`, :pull:`36`).

Internal Changes
~~~~~~~~~~~~~~~~
- Tables are removed from the package and stored in an extra github repo.
- Tables are download at first access using pooch.
- New setup structure and github ci (:pull:`16`).
- Legacy modules are ignored for coverage (:pull:`23`).
- Cleaned up import structure, avoid unneccessary warnings (:pull:`42`).

Documentation
~~~~~~~~~~~~~
- Lots of new python notebooks are rendered into the documentation.
- Documentation now includes preprocessing and pressure interpolation examples.


0.1.0 (23 July 2020)
--------------------

* First release on PyPI.
