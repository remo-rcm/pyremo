======
pyremo
======

Python tools for the regional climate model `REMO <https://www.remo-rcm.de>`_.

.. image:: https://zenodo.org/badge/282037812.svg
   :target: https://zenodo.org/badge/latestdoi/282037812

.. image:: https://github.com/remo-rcm/pyremo/actions/workflows/ci.yaml/badge.svg
    :target: https://github.com/remo-rcm/pyremo/actions/workflows/ci.yaml

.. image:: https://github.com/remo-rcm/pyremo/actions/workflows/ci-extensions.yaml/badge.svg
    :target: https://github.com/remo-rcm/pyremo/actions/workflows/ci-extensions.yaml

.. image:: https://codecov.io/gh/remo-rcm/pyremo/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/remo-rcm/pyremo

.. image:: https://img.shields.io/pypi/v/pyremo.svg
        :target: https://pypi.python.org/pypi/pyremo

.. image:: https://readthedocs.org/projects/pyremo/badge/?version=latest
        :target: https://pyremo.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://anaconda.org/conda-forge/pyremo/badges/version.svg
    :target: https://anaconda.org/conda-forge/pyremo

.. image:: https://results.pre-commit.ci/badge/github/remo-rcm/pyremo/master.svg
   :target: https://results.pre-commit.ci/latest/github/remo-rcm/pyremo/master
   :alt: pre-commit.ci status

.. image:: https://www.codefactor.io/repository/github/remo-rcm/pyremo/badge
   :target: https://www.codefactor.io/repository/github/remo-rcm/pyremo
   :alt: CodeFactor



Features
--------

* Easy access to Remo meta information and data
* API based on xarray data structures
* Pressure interpolation
* Pre and post processing
* Includes basic physics package for REMO

Installation
------------

We recommend installing ``pyremo`` with conda:

.. code-block:: console

    conda install -c conda-forge pyremo


Installation from source
------------------------

We don't recommend to pip install ``pyremo`` because some of the dependencies require pre-compiled packages
that won't work with pip. For instructions to install py-cordex from source, please have a look
at the `contributing guide <https://pyremo.readthedocs.io/en/stable/contributing.html>`_.
If you want to contribute, please get in contact as early as possible, e.g.,  using `draft pull requests <https://github.blog/2019-02-14-introducing-draft-pull-requests>`_.

Fortran extensions
------------------

There are two sub-packages that are extra private dependencies and contain Fortran extensions. For example, the preprocessing module :code:`preproc` will require the installation
of the legacy source code for preprocessing which is packaged in

* https://gitlab.dkrz.de/remo/pyintorg

For the pressure interpolation :code:`prsint`, you will need to install the additional package:

* https://gitlab.dkrz.de/remo/pydruint

Note, that you will have to install these packages from source which will require a fortran compiler (e.g. :code:`gfortran`).
If you require access to those packages, please request access to the REMO group in the DRKZ gitlab.
If you have access, you can install those extension directly from the gitlab, e.g.

.. code-block:: console

    pip install git+http://gitlab.dkrz.de/remo/pyintorg.git
    pip install git+http://gitlab.dkrz.de/remo/pydruint.git

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
