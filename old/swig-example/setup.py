from distutils.core import setup,Extension

# this one is built with swig, you'll have
# to run swig separately with
#  swig -python -c++ -o tarray_wrap.cpp tarray.i

# if you run 
#   python setup.y build_ext --inplace
# it will put it in cwd for easier testing

import numpy
include_dirs=numpy.get_include()
# can be a list
ext_modules = [Extension('_tarray', sources=['tarray_wrap.cpp',
                                             'tarray.cpp'])]
py_modules = ['tarray']


# data_files copies the ups/esutil.table into prefix/ups
setup(name='tarray',
      description='Test extension module C++ class built with swig',
      ext_modules=ext_modules,
      py_modules=py_modules,
      include_dirs=include_dirs)
