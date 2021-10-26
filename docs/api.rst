.. currentmodule:: pyremo 

#############
API reference
#############

This page provides an auto-generated summary of the pyremo API.


Top-level functions
===================

.. autosummary::
   :toctree: generated/

   remo_domain
   domain_info
   open_remo_dataset
   update_meta_info
   parse_dates


Physics
=======

.. autosummary::
   :toctree: generated/

   physics.pressure
   physics.liquid_water_content
   physics.specific_humidity
   physics.relative_humidity


Pressure interpolation
======================

.. autosummary::
   :toctree: generated/

   prsint.pressure_interpolation

   
Preprocessor
============

.. autosummary::
   :toctree: generated/

   preproc.remap
   preproc.to_netcdf
   preproc.to_tar

CF preprocessor (CMIP5/CMPI6)
-----------------------------

.. autosummary::
   :toctree: generated/

   preproc.gfile

ECMWF cmorizer (ERA5)
-------------------------

.. autosummary::
   :toctree: generated/

   preproc.ERA5
   
Converting ECMWF data
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ERA5.to_xarray
   ERA5.gfile
