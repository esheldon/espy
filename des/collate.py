"""

collate outputs

"""

from __future__ import print_function
import numpy
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

        add_dt = comb.dtype.descr + [('se1','f8'),('se2','f8'),('se','f8')]
        newdata=eu.numpy_util.add_fields(orig, add_dt)
        eu.numpy_util.copy_fields(comb, newdata)

        self.add_sellip(newdata)

        return newdata

    def add_sellip(self, data):
        data['se1']=9999
        data['se2']=9999
        data['se']=9999

        print("adding sellip")

        w,=numpy.where(data['weight'] > 0.0)

        if w.size > 0:
            print("    %d/%d with weight > 0" % (w.size,data.size))

            xxsum=data['xxsum'][w]
            xysum=data['xysum'][w]
            yysum=data['yysum'][w]

            T = xxsum + yysum
            e1 = (xxsum - yysum)/T
            e2 = 2*xysum/T
            e=numpy.sqrt(e1**2 + e2**2)

            data['se1'][w] = e1
            data['se2'][w] = e2
            data['se'][w] = e
