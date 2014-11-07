"""

collate outputs

"""

from __future__ import print_function
from .files import *

class Collator(dict):
    def __init__(self, run):
        conf=cascade_config(run)
        self.update(conf)

    def collate(self):
        """
        collate the combined file with the original catalog
        """
        import fitsio

        d=get_collated_dir(self['run'])
        fname=get_collated_file(self['run'])

        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        data=self.get_collated_data()

        print("writing:",fname)
        with fitsio.FITS(fname,'rw',clobber=True) as fits:
            fits.write(data)
        
    def get_collated_data(self):
        """
        read orig and xshear output and combine them
        """
        import esutil as eu 

        orig=read_lcat_original(self['lens_conf']['lcat_name'])
        comb=read_combined(self['run'])

        print("matching")
        index_col=self['lens_conf']['index_col']
        mo,mc=eu.numpy_util.match(orig[index_col], comb['index'])

        print("    matched %d/%d" % (mo.size, orig.size))

        orig=orig[mo]
        comb=comb[mc]

        print("collating")
        newdata=eu.numpy_util.add_fields(orig, comb.dtype.descr)
        eu.numpy_util.copy_fields(comb, newdata)

        return newdata
