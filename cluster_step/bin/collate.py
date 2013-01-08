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

parser.add_option('-r','--run',default=None,
                  help='The run id, required')

parser.add_option('-f','--field',default='Ts2n',
                  help="field for S/N, default %default")

parser.add_option('--s2n',default=20, help=("threshold in s/n"))

parser.add_option('--s2',default=1.0,
                  help='restrict s2 less than this value, default %d')


class Collator(object):
    def __init__(self, run, **keys):
        self.run=run

        self.s2_max=keys['s2_max']
        self.s2n_min=keys['s2n_min']
        self.s2n_field=keys['s2n_field']
        self.objtype='gexp'

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
                                          shnum=shnum,
                                          s2_max=self.s2_max,
                                          s2n_min=self.s2n_min)
        dir=os.path.dirname(path)
        if not os.path.exists(dir):
            print 'making directory:',dir
            os.makedirs(dir)

        print 'writing:',path
        with fitsio.FITS(path,mode='rw',clobber=True) as fits:
            fits.write(exp_data)

    def set_use(self, data):
        logic1 = (data['s2'] < self.s2_max)
        logic1 = logic1 & (data[self.s2n_field] > self.s2n_min)

        logic2=logic1 & (data['objtype'] == self.objtype)

        w1,=where(logic1)
        w2,=where(logic2)

        frac1=float(w1.size)/data.size
        frac2=float(w2.size)/data.size
        print 'keep1 %s/%s %.2f' % (w1.size,data.size,frac1)
        print 'keep2 %s/%s %.2f' % (w2.size,data.size,frac2)

        if w1.size > 0:
            data['use1'][w1]=1
        if w2.size > 0:
            data['use2'][w2]=1


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
            ('use2','i2')]

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

    if options.run is None:
        parser.print_help()
        sys.exit(1)

    run=options.run
    s2_max=float(options.s2)
    s2n_min=float(options.s2n)
    s2n_field=options.field

    collator=Collator(run, 
                      s2_max=s2_max,
                      s2n_min=s2n_min,
                      s2n_field=s2n_field)
    collator.go()

main()
