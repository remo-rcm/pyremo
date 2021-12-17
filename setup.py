#!/usr/bin/env python
from setuptools import setup, find_packages

# get version
with open("pyremo/version.py") as f:
    line = f.readline().strip().replace(" ", "").replace('"', "")
    version = line.split("=")[1]
    __version__ = version

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read().strip().split("\n")

setup(
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Common REMO python tools",
    entry_points={
        'console_scripts': [
            'prsint=pyremo.cli:main_prsint',
            'remo_analysis=pyremo.cli:main_analysis'
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyremo',
    name='pyremo',
    packages=find_packages(include=['pyremo',
                                    'pyremo.analysis']),
    test_suite='tests',
    url='https://github.com/remo-rcm/pyremo',
    version=__version__,
    zip_safe=False,
)
