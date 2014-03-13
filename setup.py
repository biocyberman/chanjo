#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from setuptools import setup, find_packages
import chanjo

# Shortcut for publishing to Pypi
# Source: https://github.com/kennethreitz/tablib/blob/develop/setup.py
if sys.argv[-1] == "publish":
  os.system("python setup.py sdist upload")
  sys.exit()

with open('README.rst') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

with open('requirements.txt') as f:
  requirements = [line.rstrip() for line in f.readlines()]

with open('dev_requirements.txt') as f:
  dev_requirements = [line.rstrip() for line in f.readlines()]

setup(
  # Metadate for upload to Pypi
  name=chanjo.__name__,
  version=chanjo.__version__,
  description=chanjo.__description__,
  long_description=readme,
  keywords='coverage sequencing clinical exome',
  platform='UNIX',
  author=chanjo.__author__,
  author_email=chanjo.__email__,
  url=chanjo.__url__,
  download_url='https://github.com/robinandeer/chanjo',
  license=license,
  packages=find_packages(exclude=('tests', 'docs')),
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Programming Language :: Python :: 3.3',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ],

  # Executable command line utilities
  scripts=[
    'chanjo.py'
  ],

  # Runtime dependencies
  install_requires = requirements,

  # Testing dependencies
  tests_require = [
    'pytest'
  ]
)
