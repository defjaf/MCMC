from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np                           # <---- New line

ext_modules = [Extension("blackbody", ["blackbody.pyx"])]

setup(
  name = 'Blackbody module',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include()],         # <---- New line
  ext_modules = ext_modules
)


### easiest way to use this is to run python setup.py build; mv --force build/lib*/blackbody.so ../submmSED/
