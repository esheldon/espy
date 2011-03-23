from __future__ import print_function
import numpy
import sdsspy
import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1
import os
import columns
import datetime

from . import stomp_maps

class Collator:
    """

    Collate the regauss outputs into columns.

    Note we put both types ('gal','star') into the same run.  This means if you
    make indices you'll have to re-make them when the other type is added.

    """
    
    def __init__(self, type='gal'):
        if type != 'gal':
            raise ValueError("add support for 'star' type")

        self.type=type
        self.minscore=0.1
        self.flags = sdsspy.flags.Flags()
        self.masktypes = ['basic','good','tycho']

    def run(self, runs=None):
        """
        Collate the collated files from the run into a columns directory.
        Multi-band columns are split to, e.g., flux_r, flux_i, ...
        """

        c = self.open_columns()
        print("Will write in coldir:", c.dir)
        if c.dir_exists():
            raise ValueError("Columns already exist, start from scratch")
        c.create()


        # create the meta column from input data
        meta = self.make_meta()
        c.write_column('meta', meta)

        if runs is None:
            win = sdsspy.window.Window()
            runs, reruns = win.runlist(self.minscore)

        nrun=len(runs) 
        print("Processing",nrun,"runs\n")

        for i in xrange(nrun):
            run=runs[i]
            print("\nProcessing run: %s  %s/%s" % (run,i+1,nrun))
            print("-"*70)
            for camcol in [1,2,3,4,5,6]:
                print("  %06i-%i" % (run,camcol))
                tmp = sdsspy.read('calibobj.'+self.type, run, camcol, 
                                  lower=True, ensure_native=True)

                if len(tmp) == 0:
                    raise ValueError("sweep is empty")
                else:
                    print("    Found",tmp.size)
                    print("    Getting primary")
                    w=self.get_primary(tmp)
                    print("    Found",w.size)
                    if w.size > 0:
                        tmp = tmp[w]
                        print("      Getting output")
                        out_dict = self.create_output(tmp)
                        print("      Writing columns")
                        for name in out_dict:
                            c.write_column(name, out_dict[name])
        print("Done")

    def get_primary(self, objs):
        primary = self.flags.val('resolve_status','survey_primary')

        w=where1((objs['resolve_status'] & primary) != 0)
        return w


    def create_output(self, st):
        bands = ['u','g','r','i','z']

        out_dict = {}
        for d in st.dtype.descr:
            name = str( d[0] )

            if len(d) == 3:
                if d[2] == 5:
                    # ignoring the higher dim stuff
                    for bandi in xrange(5):
                        fname = name+'_'+bands[bandi]
                        out_dict[fname] = st[name][:, bandi]
            else:
                out_dict[name] = st[name]

        cmodelmag_dered = sdsspy.make_cmodelmag(st, doerr=False, dered=True)
        modelmag = sdsspy.nmgy2mag(st['modelflux'])
        modelmag_dered = modelmag - st['extinction']

        for f in sdsspy.FILTERCHARS:
            fnum = sdsspy.FILTERNUM[f]
            out_dict['cmodelmag_dered_'+f] = cmodelmag_dered[:,fnum]
            out_dict['modelmag_dered_'+f] = modelmag_dered[:,fnum]

        out_dict['photoid'] = sdsspy.photoid(st)

        self.set_maskflags(out_dict)
        return out_dict

    def create_indices(self):
        c = self.open_columns()

        cnames = ['cmodelmag_dered_r','modelmag_dered_r',
                  'run','camcol','ra','dec','thing_id','photoid',
                  'inbasic','ingood','intycho']

        for n in cnames:
            print("Creating index for:",n)
            c[n].create_index(force=True)



    def make_meta(self):
        meta={}
        meta['photo_sweep'] = os.environ['PHOTO_SWEEP']

        now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
        meta['date'] = now
        return meta

    def open_columns(self):
        dir = self.columns_dir()
        return columns.Columns(dir)

    def columns_dir(self):
        self.photo_sweep=os.environ['PHOTO_SWEEP']
        dirname=os.path.basename(self.photo_sweep)
        dirname += '_new'
        dir=os.environ['SWEEP_REDUCE']
        dir = path_join(dir,dirname)
        if not os.path.exists(dir):
            os.makedirs(dir)
                        
        dir=path_join(dir ,'prim'+self.type+'.cols')
        return dir


    
    def set_maskflags(self, new):
        self.load_masks()
        ra = new['ra']
        dec = new['dec']
        print("        Checking mask: ",end='')
        inside = numpy.zeros(ra.size, dtype='i1')
        for masktype in self.masktypes:
            name = 'in'+masktype
            print(" ",name,end='')

            inside[:] = self.masks[masktype].Contains(ra,dec,'eq')
            new[name] = inside.copy()
     
        print("")



    def load_masks(self):
        if not hasattr(self,'masks'):
            self.masks = {}
            for masktype in self.masktypes:
                self.masks[masktype] = stomp_maps.load('boss',masktype)




