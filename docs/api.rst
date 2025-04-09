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

Preconfigures preprocessors
---------------------------

.. autosummary::
   :toctree: generated/

   preproc.Preprocessor
   preproc.ERA5Preprocessor
   preproc.CFPreprocessor
   preproc.RemoPreprocessor

Preprocessing functions
-----------------------

.. autosummary::
   :toctree: generated/

   preproc.get_gcm_dataset
   preproc.get_gcm_gfile
   preproc.remap
   preproc.remap_remo
   preproc.ERA5
   preproc.ERA5.gfile

Tutorial
========

.. autosummary::
   :toctree: generated/

   tutorial.open_dataset
   tutorial.load_dataset
