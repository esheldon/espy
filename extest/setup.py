#from distutils.core import setup, Extension
#setup(name="cosmolib", version="1.0",
#      ext_modules=[Extension("cosmolib", ["cosmolib_pywrap.c","cosmolib.c"])])


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    #config = Configuration('cosmology',parent_package,top_path)
    config = Configuration(None,parent_package,top_path)
    sources=["cosmolib_pywrap.c","cosmolib.c"]
    config.add_extension('cosmolib', sources)

    return config

if __name__=="__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
