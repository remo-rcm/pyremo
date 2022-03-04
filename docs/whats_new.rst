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

 
Internal Changes
~~~~~~~~~~~~~~~~
- Added fortran extensions ``pyintorg`` and ``pydruint`` to extra dependencies.
- Added ``CI-extensions.yaml`` for git actions testing with fortran extensions.

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
