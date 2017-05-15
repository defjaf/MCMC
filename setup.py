# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import re

## from https://gehrcke.de/2014/02/distributing-a-python-command-line-application/
## note that it requires double quotes around the version number... should fix!
version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('MCMC/__init__.py').read(),
    re.M
    ).group(1)

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

setup(name='MCMC',
      version=version,
      description='MCMC driver.',
      url='http://github.com/defjaf/MCMC',
      author='Andrew H. Jaffe',
      author_email='a.h.jaffe@gmail.com',
      install_requires = ['numpy', 'matplotlib', 'scipy', 'progressbar-ipython'],
      packages=find_packages(exclude=('test',)),
      package_data={'MCMC': [
        'binnedCl/*.txt', 
        'ClData/*.dat', 
      ##  'topology/likelihood/*/dat/*.dat',    ## don't include for now, j
                                                ##  just leave in source dir
        'topology/likelihood/*/*.f',
        'topology/likelihood/*/*.for',
        'topology/likelihood/wmap/*',
        'submmSED/dat/*',
        'submmSED/M31/*'
    ]}
)

####  testing reminder 
# virtualenv --python=python2 venvpy27
# source venvpy27/bin/activate
# python setup.py install

# virtualenv --python=python3 venvpy3
# source venvpy3/bin/activate
# python setup.py install
