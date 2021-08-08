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
  
Internal Changes
~~~~~~~~~~~~~~~~
- Tables are removed from the package and stored in an extra github repo.
- Tables are download at first access using pooch.
