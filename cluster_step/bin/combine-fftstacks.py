"""
    %prog [options]
"""

import sys
import os
import numpy
import lensing
from numpy import zeros, sqrt, where

import cluster_step
from cluster_step import files
from cluster_step import sh1exp, sh2exp

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-s','--shnum',default=None,
                  help='The shear number, required')

parser.add_option('-p','--psfnum',default=None,
                  help='The psf number, required')

class Combiner(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if (options.run is None 
                or options.psfnum is None
                or options.shnum is None):
            parser.print_help()
            sys.exit(1)

        self.run=options.run
        self.psfnum=int(options.psfnum)
        self.shnum=int(options.shnum)

    def go(self):
        import fitsio
        import esutil as eu

        reslist=[]

        nstack=0
        npsf_stack=0
        psf_skyvar=0.0
        for ccd in xrange(1,62+1):
            path=files.get_output_path(ftype='shear',
                                       run=self.run,
                                       psfnum=self.psfnum,
                                       shnum=self.shnum,
                                       ccd=ccd)
            
            print path
            with fitsio.FITS(path) as fits:

                res=fits[1].read()

                psf_stacks_real=fits[2].read()
                psf_stacks_imag=fits[3].read()

                psf_h=fits[2].read_header()

                npsf_stack += psf_h['nstack']
                psf_skyvar += psf_h['skyvar']

                im_stacks=fits[4].read()

            reslist.append(res)
            if ccd==1:
                psf_stack_real_tot = psf_stacks_real.copy()
                psf_stack_imag_tot = psf_stacks_imag.copy()
                im_stacks_tot = im_stacks.copy()
            else:
                psf_stack_real_tot += psf_stacks_real
                psf_stack_imag_tot += psf_stacks_imag
                # these sum across multiple dimensions
                im_stacks_tot['images_real'] += im_stacks['images_real']
                im_stacks_tot['images_imag'] += im_stacks['images_imag']
                im_stacks_tot['nstack'] += im_stacks['nstack']
                im_stacks_tot['skyvar'] += im_stacks['skyvar']
        
        res=eu.numpy_util.combine_arrlist(reslist)

        self._write_data(res, 
                         psf_stack_real_tot, psf_stack_imag_tot, npsf_stack, psf_skyvar,
                         im_stacks_tot, psf_h['boost'])

    def _write_data(self, res, psf_stack_real, psf_stack_imag, npsf_stack, psf_skyvar,
                    im_stacks, boost):
        import fitsio
        path=files.get_output_path(ftype='shear-stack',
                                   run=self.run,
                                   psfnum=self.psfnum,
                                   shnum=self.shnum)


        psf_header={'skyvar':psf_skyvar,
                    'nstack':npsf_stack,
                    'boost':boost}
        print path
        with fitsio.FITS(path, mode='rw', clobber=True) as fits:
            fits.write(res)
            fits.write(psf_stack_real, header=psf_header)
            fits.write(psf_stack_imag, header=psf_header)
            fits.write(im_stacks)
    
 

def main():
    c=Combiner()
    c.go()

main()
