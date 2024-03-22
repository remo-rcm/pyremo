.. currentmodule:: pyremo

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo

v0.7.0 (22 March 2024)
----------------------

This release makes the codes API public. It also completely retires the cmor module in favour of `py-cordex.cmor <https://py-cordex.readthedocs.io/en/latest/generated/cordex.cmor.cmorize_variable.html>`_.

New Features
~~~~~~~~~~~~

- Added :py:meth:`codes.search` for searching in the code table using keyword arguments (:pull:`160`).

Internal Changes
~~~~~~~~~~~~~~~~

- Pin ``numpy < 1.26.3``, causes problems with installing Fortran extensions (:pull:`158`).
- Update ERA5 cmorizer and preprocessor (:pull:`146`).

Breaking Changes
~~~~~~~~~~~~~~~~

- CMOR module is retired and removed in favour of `py-cordex <https://py-cordex.readthedocs.io/en/latest/generated/cordex.cmor.cmorize_variable.html#cordex-cmor-cmorize-variable>`_ (:pull:`157`).

v0.6.1 (20 October 2023)
------------------------

Patch release to fix preprocessing issues.

Internal Changes
~~~~~~~~~~~~~~~~

- Fix SST interpolation and extrapolation for CF preprocessor (:pull:`140`).
- Fix ERA5 cmorizer subprocess call (:pull:`141`).

v0.6.0 (25 August 2023)
-----------------------

This release contains some internal refactoring for upcoming CMIP6 simulations. There has also been an update on the
documentation. The :py:meth:`preproc.ERA5` version is now compatible with the latest conventions of the
`DKRZ ECMWF data pool <https://docs.dkrz.de/doc/dataservices/finding_and_accessing_data/era_data/>`_. There
has also been some significant clean up of the internal preprocessor modules.

Internal Changes
~~~~~~~~~~~~~~~~

- Fixed E721 rule, switched to `ruff <https://github.com/astral-sh/ruff>`_ for linting in pre-commit, added ``dependabot.yml`` (:pull:`132`).
- Renamed `output_pattern` to :py:meth:`file_pattern` and added afile patterns (:pull:`132`).
- Domain notebook running on readthedocs (:pull:`126`).
- Refactoring of core modules, now using `domain_id` keyword instead of `short_name` (:pull:`125`).
- Renamed ``master`` to ``main``.

Deprecations
~~~~~~~~~~~~

- `output_pattern` is deprecated, please use to :py:meth:`file_pattern` instead.

Breaking Changes
~~~~~~~~~~~~~~~~

- Refactoring of :py:meth:`preproc.ERA5` for new DKRZ data pool conventions (:pull:`129`).
- Refactoring and clean up of preprocessing module (:pull:`130`).

v0.5.1 (24 May 2023)
--------------------

Updates for REMO run workflows.

Bugfixes
~~~~~~~~

- Fixed pattern for f-files and p-files (:pull:`117`, :pull:`122`).
- Fixed pressure level sorting (:pull:`120`).
- Added `z_coord` keyword to :py:meth:`physics.pressure` (:pull:`121`).

v0.5.0 (3 March 2023)
---------------------

This release comes with updates for CMIP6 and ERA5 downscaling.

New Features
~~~~~~~~~~~~

- Added :py:meth:`output_pattern` for creating output file naming patterns (:pull:`113`).

Internal Changes
~~~~~~~~~~~~~~~~

- Updates for CI pipeline (:pull:`108`).
- Updates for ERA5 cmorizer to work with the new DKRZ ERA5 catalog (:pull:`112`).
- UV correction is now optional in :py:meth:`cmor.remap_remo` (:pull:`111`).

Documentation
~~~~~~~~~~~~~

- Updates GHG notebook with ssp scenario table (:pull:`114`).


v0.4.1 (19 January 2023)
------------------------

Patch release for preprocessor updates.

Internal Changes
~~~~~~~~~~~~~~~~

- Allow multiple input directories for CF preprocessor (:pull:`106`).
- Added ``scratch`` keyword to gfile API (:pull:`106`).
- Drop ``python3.7`` support in favour of latest xarray releases (:pull:`106`).


v0.4.0 (15 November 2022)
-------------------------

This release comes with an updated preprocessing module for preprocessing of CMIP6 model data and the double nesting preprocessor (includes `ptop`). There are also some additional tools for working with forcing files and an updated documentation.


New Features
~~~~~~~~~~~~

- Added :py:meth:`magic_numbers` to API (:pull:`85`).
- New module for gfile creation from CMIP6 datasets (:pull:`78`, :pull:`89`).
- Additional CF preprocessor option which ueses `xesmf` for horizontal interpolation (:pull:`100`). This is not part of the user API yet.
- New notebook on greenhouse gas concentration (GHG). This documents how we create GHG forcing tables for CMIP6 downscaling (:pull:`92`).
- Command line interface ``pradd-vars`` for variable replacement: This tool can be used to replace soil variables in a forcing file with data from a REMO output file (*warm soil*) (:pull:`93`).
- Updates for (double-nesting) preprocessor (:pull:`98`), includes an implementation for using a top level pressure (`ptop`) during vertical interpolation (:pull:`104`).

Internal Changes
~~~~~~~~~~~~~~~~

- Pinned ``setuptools < 60.0`` due to Fortran build system (:pull:`94`).
- Table fetching now ignores file hashes. The tables will now be pulled from the ``main`` branch instead of ``master`` in the `tables repository <https://github.com/remo-rcm/tables>`_ (:pull:`97`). Older versions will rely on the ``master`` branch, so that one is now frozend and protected.
- Updates and fixes for the preprocessing module (:pull:`100`, :pull:`104`).

Documentation
~~~~~~~~~~~~~

- Documentation updates (:pull:`87`).


v0.3.4 (11 July 2022)
---------------------

Patch release to use ``pyremo`` for cmorization tests with `CORDEX-CMIP6 cmor tables <https://github.com/WCRP-CORDEX/cordex-cmip6-cmor-tables>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Updated ``pyremo.cmor`` API to choose cmor grids table (:pull:`83`)

v0.3.3 (27 June 2022)
---------------------

Patch release to `fix version bug <https://github.com/remo-rcm/pyremo/commit/0ab457d7ba6f828497059797f66d218d26ca954a>`_.

v0.3.2 (24 June 2022)
---------------------

Patch release to fix preprocessor bugs.

Internal Changes
~~~~~~~~~~~~~~~~

- Switched to automatic version numbering using ``setuptools_scm``, added ``publish-pypi.yaml`` workflow (:pull:`80`).


Bugfixes
~~~~~~~~

- Updated path resource for ``pyremo.analysis`` for use at DKRZ on levante filesystem (:pull:`76`).
- Fixed SST interpolation in preprocessing of CMIP models, uses now masks with xesmf (:pull:`79`).
- Fixed global attributes in forcing data inherited from driving model (:pull:`81`).


v0.3.1 (4 May 2022)
-------------------

Patch release to fix ``.zenodo.json``, update documentation and add pre-commit hooks.

Internal Changes
~~~~~~~~~~~~~~~~

- Added ``nbQA`` to pre-commit hooks and formatted notebooks (:pull:`72`).

Bugfixes
~~~~~~~~

- Fixed ``.zenodo.json`` format and added pre-commit hook (:pull:`71`).


v0.3.0 (4 May 2022)
-------------------

This release comes with a new `double nesting preprocessing API <https://pyremo.readthedocs.io/en/latest/generated/pyremo.preproc.remap_remo.html#pyremo.preproc.remap_remo>`_,
a command line interface for the `pressure interpolation <https://pyremo.readthedocs.io/en/latest/prsint.html>`_, updated
`cmorization <https://pyremo.readthedocs.io/en/latest/cmorization.html>`_ features and a new ERA5 cmorizer to handle ERA5.1 data.
Documentation has also improved and contains an updated tutorial on how to `explore REMO NetCDF output <https://pyremo.readthedocs.io/en/latest/remo-dataset.html>`_.

New Features
~~~~~~~~~~~~

- Added double nesting preprocessor API (:pull:`34`).
- ``prsint``: a command line interface for the pressure interpolator (:pull:`53`, :pull:`66`, :pull:`68`).
- ERA5 cmorizer includes ERA5.1 to work with DKRZ data pool (:pull:`52`).
- Added ``python3.10`` support (:pull:`46`).
- Updated cmorization module for use with CMIP6 (:pull:`48`, :pull:`49`, :pull:`51`) and CORDEX vocabulary. The underlying tables are only used for testing and should not yet be uses for actual data publication, see also `here <https://github.com/euro-cordex/py-cordex/pull/55>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Added tests for `pyremo.preproc` module (:pull:`34`).
- Updated documentation including an improved contribution guide (:pull:`63`, :pull:`65`).
- Fixed issues with ``dask.delayed`` in ``preproc.era5`` running in batch mode (:pull:`58`).
- Added ``.pre-commit-config.yaml`` and ``linting.yaml`` to run linter checks (:pull:`55`).
- Code is reformatted to apply to ``flake8`` and ``black`` conventions (:pull:`55`).
- ERA5 cmorizer works with pandas datatables (:pull:`52`).
- Tutorial `data source <https://github.com/remo-rcm/pyremo-data>`_ is now on github (:pull:`47`).
- Added fortran extensions ``pyintorg`` and ``pydruint`` to extra dependencies (:pull:`64`).
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
