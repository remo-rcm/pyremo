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
   parse_absolute_time
   preprocess
   magic_number
   magic_numbers
   file_pattern

Codes
=====

.. autosummary::
   :toctree: generated/

   codes.get_dict
   codes.search

Physics
=======

.. autosummary::
   :toctree: generated/

   physics.pressure
   physics.liquid_water_content
   physics.specific_humidity
   physics.relative_humidity
   physics.precipitation_flux
   physics.water_vapour
   physics.specific_humidity_from_dewpoint
   physics.relative_humidity_from_dewpoint
   physics.surface_runoff_flux
   physics.surface_downwelling_shortwave_flux_in_air
   physics.surface_downwelling_longwave_flux_in_air
   physics.toa_incoming_shortwave_flux
   physics.water_evapotranspiration_flux


Pressure interpolation
======================

.. autosummary::
   :toctree: generated/

   prsint.pressure_interpolation


REMO preprocessor
=================

CMIP (CF) preprocessor
----------------------

.. autosummary::
   :toctree: generated/

   preproc.gfile
   preproc.remap
   preproc.to_netcdf
   preproc.to_tar

Double nesting preprocessor
---------------------------

.. autosummary::
   :toctree: generated/

   preproc.remap_remo


ECMWF cmorizer (ERA5)
---------------------

.. autosummary::
   :toctree: generated/

   preproc.ERA5


Converting ECMWF data
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   preproc.ERA5.gfile


Tutorial
========

.. autosummary::
   :toctree: generated/

   tutorial.open_dataset
   tutorial.load_dataset
