from distutils.core import setup, Extension

import numpy
include_dirs=[numpy.get_include()]
setup(name="_cosmolib", version="1.0",
      ext_modules=[Extension("_cosmolib", ["cosmolib_pywrap.c","cosmolib.c"])],
      include_dirs=include_dirs)


"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    #config = Configuration('cosmology',parent_package,top_path)
    config = Configuration(None,parent_package,top_path)
    sources=["cosmolib_pywrap.c","cosmolib.c"]
    config.add_extension('_cosmolib', sources)

    return config

if __name__=="__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
"""
