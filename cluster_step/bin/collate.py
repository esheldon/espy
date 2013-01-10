"""
    %prog [options]

The use1 field is set based on s2n_min and s2_max.

The use2 is an and of use1 and objtype=='gexp'
"""

import sys
import os
from numpy import where, zeros, sqrt
import cluster_step
from cluster_step import files, stats
from esutil.numpy_util import match, strmatch, combine_arrlist, copy_fields

import fitsio

from optparse import OptionParser

parser=OptionParser(__doc__)

class Collator(object):
    def __init__(self, run, **keys):
        self.run=run
        self.conf=files.read_config(self.run)

    def go(self):
        for psfnum in files.PSFNUMS:
            for shnum in files.SHNUMS:
                
                exp_data=[]
                for ccd in files.CCDS:
                    pipe=self.load_single(psfnum,shnum,ccd)

                    ccd_data_i=self.get_matched_struct(pipe)

                    self.set_use(ccd_data_i)

                    exp_data.append( ccd_data_i )

                exp_data=combine_arrlist(exp_data)

                self.write_data(exp_data, psfnum, shnum)

    def write_data(self, exp_data, psfnum, shnum):
        path=files.get_julia_collate_path(run=self.run,
                                          psfnum=psfnum,
                                          shnum=shnum)
        dir=os.path.dirname(path)
        if not os.path.exists(dir):
            print 'making directory:',dir
            os.makedirs(dir)

        print 'writing:',path
        with fitsio.FITS(path,mode='rw',clobber=True) as fits:
            fits.write(exp_data)

    def set_use(self, data):
        sratio=sqrt(1./data['s2'])

        isexp=(data['objtype'] == 'gexp')
        isdev=(data['objtype'] == 'gdev')
        
        #size1 = (sratio > 1.0)
        #size2 = (sratio > sqrt(2))
        Ts2n_20 = (data['Ts2n'] > 20)

        exp_or_dev20 = ( isexp | (isdev & Ts2n_20) )
        exp20 = isexp & Ts2n_20
        
        # objsize > psfsize and either exp or well measured dev
        logic1 = exp_or_dev20

        # objsize > psfsize and is exp
        logic2 = isexp

        # objsize > psfsize and and well measured
        logic3 = Ts2n_20

        # objsize > psfsize and and well measured
        logic4 = isexp & Ts2n_20


        w1,=where(logic1)
        w2,=where(logic2)
        w3,=where(logic3)
        w4,=where(logic4)

        frac1=float(w1.size)/data.size
        frac2=float(w2.size)/data.size
        frac3=float(w3.size)/data.size
        frac4=float(w4.size)/data.size
        print 'keep1 %s/%s %.2f' % (w1.size,data.size,frac1)
        print 'keep2 %s/%s %.2f' % (w2.size,data.size,frac2)
        print 'keep3 %s/%s %.2f' % (w3.size,data.size,frac3)
        print 'keep4 %s/%s %.2f' % (w4.size,data.size,frac4)

        if w1.size > 0:
            data['use1'][w1]=1
        if w2.size > 0:
            data['use2'][w2]=1
        if w3.size > 0:
            data['use3'][w3]=1
        if w4.size > 0:
            data['use4'][w4]=1


    def get_matched_struct(self, pipe):
        #dt=pipe.cat.dtype.descr
        #dt += [
        dt=[
            ('simid','i4'),
            ('ccd','i2'),
            ('id','i4'),
            ('row','f8'),
            ('col','f8'),
            ('s2n','f8'),
            ('Ts2n','f8'),
            ('s2','f8'),
            ('sratio','f8'),
            ('objtype','S4'),
            ('g','f8',2),
            ('gcov','f8',(2,2)),
            ('gsens','f8',2),
            ('weight','f8'),
            ('use1','i2'),
            ('use2','i2'),
            ('use3','i2'),
            ('use4','i2')]

        # note this s2n min is just the limit used in the
        # shear measurement code
        wgal=pipe.get_gals(s2n_min=self.conf['shear_s2n_min'])

        gals=pipe.cat[wgal]
        output=zeros(gals.size, dtype=dt)

        print wgal.size, pipe.shear_res.size
        w,=where(gals['simid'] != pipe.shear_res['simid'])
        if w.size != 0:
            raise ValueError("gals and shear don't line up")

        mcat,mshear=match(gals['simid'], pipe.shear_res['simid'])
        if mshear.size != pipe.shear_res.size:
            mess="""not all shear objects matched
                by simid: %d %d""" % (pipe.shear_res.size,mshear.size)
            print mess
         
        # gets simid, row, col
        copy_fields(gals, output)

        output['ccd'] = pipe['ccd']

        output['s2n'][:] = pipe.shear_res['s2n_w'][:]
        output['Ts2n'][:] = pipe.shear_res['Ts2n'][:]
        output['s2'][:] = pipe.shear_res['s2'][:]
        output['sratio'][:] = sqrt(1./output['s2'][:])

        for i in xrange(output.size):
            output['objtype'][i] = pipe.shear_res['model'][i]

        output['g'][:,:] = pipe.shear_res['g'][:,:]
        output['gcov'][:,:,:] = pipe.shear_res['gcov'][:,:,:]
        output['gsens'][:,:] = pipe.shear_res['gsens'][:,:]

        wts=stats.get_weights(pipe.shear_res['gcov'][:,:,:])
        output['weight'][:] = wts[:]

        return output


    def load_single(self,psfnum,shnum,ccd):
        pipe=cluster_step.pipe.Pipe(run=self.run, psfnum=psfnum, shnum=shnum,
                                    ccd=ccd)
        if not hasattr(pipe,'shear_res'):
            mess="""
missing:
    %s
    %s
    %s
    %s
            """ % (self.run,psfnum,shnum,ccd)
            raise ValueError(mess)

        return pipe

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run=args[0]

    collator=Collator(run)
    collator.go()

main()
