.. currentmodule:: pyremo

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo 

v0.2.0 (unreleased)
-------------------

This is a major restructuring release. 

New Features
~~~~~~~~~~~~
- Included physics package with xarray API.
- Included preprocessing interface for xarray data structures.
- Included ERA5 cmorizer.
- Includes experimental cmorization module (:pull:`20`, :pull:`21`, :pull:`22`,
  :pull:`23`, :pull:`25`, :pull:`26`).
- Included production analysis (:pull:`33`).
 
Internal Changes
~~~~~~~~~~~~~~~~
- Tables are removed from the package and stored in an extra github repo.
- Tables are download at first access using pooch.
- New setup structure and github ci (:pull:`16`).
- Legacy modules are ignored for coverage (:pull:`23`).

Documentation
~~~~~~~~~~~~~
- Lots of new python notebooks are rendered into the documentation. 
- Documentation now includes preprocessing and pressure interpolation examples.


0.1.0 (2020-07-23)
------------------

* First release on PyPI.
