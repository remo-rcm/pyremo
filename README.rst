======
pyremo
======

.. image:: https://zenodo.org/badge/282037812.svg
   :target: https://zenodo.org/badge/latestdoi/282037812

.. image:: https://github.com/remo-rcm/pyremo/actions/workflows/ci.yaml/badge.svg
    :target: https://github.com/remo-rcm/pyremo/actions/workflows/ci.yaml

.. image:: https://codecov.io/gh/remo-rcm/pyremo/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/remo-rcm/pyremo

.. image:: https://img.shields.io/pypi/v/pyremo.svg
        :target: https://pypi.python.org/pypi/pyremo

.. image:: https://readthedocs.org/projects/pyremo/badge/?version=latest
        :target: https://pyremo.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status
        
.. image:: https://anaconda.org/conda-forge/pyremo/badges/installer/conda.svg
    :target: https://anaconda.org/conda-forge/pyremo

.. image:: https://www.codefactor.io/repository/github/remo-rcm/pyremo/badge
   :target: https://www.codefactor.io/repository/github/remo-rcm/pyremo
   :alt: CodeFactor



Common Remo python tools

* Free software: MIT license
* Documentation: https://pyremo.readthedocs.io.

Features
--------

* Easily access to Remo meta information and data
* API based on xarray data structures
* Pressure interpolation
* Pre and post processing
* Includes basic physics package

Installation
------------

You can install the package directly from github using pip:


.. code-block:: console

    pip install git+https://github.com/remo-rcm/pyremo


If you want to contribute, I recommend cloning the repository and installing the package in development mode, e.g.


.. code-block:: console

    git clone https://github.com/remo-rcm/pyremo
    cd pyremo
    pip install -e .


This will install the package but you can still edit it and you don't need the package in your :code:`PYTHONPATH`

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
