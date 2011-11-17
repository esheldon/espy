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

def columns_dir(type):
    collator = Collator(type)
    return collator.columns_dir()
def open_columns(type):
    collator = Collator(type)
    print("opening columns:",collator.columns_dir())
    return collator.open_columns()

def match_columns(photoid, tags, type='primgal'):
    """

    Take the input photoids and match it to the input photo columns.  Extract
    the requested columns.

    """
    if isinstance(tags, str):
        tags=[tags]
    cols = open_columns(type)
    add_descr = []
    for tag in tags:
        if tag not in cols:
            raise ValueError("Requested tag not in database:",tag)

        print("will extract:",cols[tag].dtype.descr[0])
        add_descr.append(cols[tag].dtype.descr[0])

    print("reading photoid from columns")
    pid = cols['photoid'][:]
    
    print("matching")
    minput, mcols = eu.numpy_util.match(photoid, pid)

    if minput.size != photoid.size:
        raise ValueError("Not all objects matched: %d/%d" % (minput.size, photoid.size))

    print("verifying")
    tpid = cols['photoid'][mcols]
    wbad=where1(tpid != photoid[minput])
    if wbad.size != 0:
        raise ValueError("extracted photoid did not match up")

    struct = cols.read_columns(tags, rows=mcols)
    return struct

class Collator:
    """

    Collate the regauss outputs into columns.

    Note we put both types ('gal','star') into the same run.  This means if you
    make indices you'll have to re-make them when the other type is added.

    """
    
    def __init__(self, type='gal'):
        #if type not in ['gal','primgal']:
        #    raise ValueError("add support for 'star' type")

        self.type=type
        if type in ['gal','primgal']:
            self.sweep_type = 'gal'
        elif type in ['star','primstar']:
            self.sweep_type = 'star'
        else:
            raise ValueError("type should be 'gal','primgal','star','primstar'")

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
            raise ValueError("coldir already exist, start from scratch")
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
                tmp = sdsspy.read('calibobj.'+self.sweep_type, run, camcol, 
                                  lower=True, ensure_native=True)

                if len(tmp) == 0:
                    raise ValueError("sweep is empty")

                print("    Found",tmp.size)
                if self.type in ['primgal','primstar']:
                    print("    Selecting survey_primary")
                    w=self.get_primary_indices(tmp)
                else:
                    print("    Selecting run_primary")
                    w=self.get_primary_indices(tmp, run_primary=True)

                print("    Found",w.size)
                if w.size > 0:
                    tmp = tmp[w]

                    print("      Getting output")
                    out_dict = self.create_output(tmp)
                    print("      Writing columns")
                    for name in out_dict:
                        c.write_column(name, out_dict[name])
        print("Done")


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

        if 'devflux' in st.dtype.names:
            #cmodelflux, cmodelflux_ivar = sdsspy.make_cmodelflux(st, doivar=True)
            cmodelmag_dered, cmodelmag_dered_err = sdsspy.make_cmodelmag(st, dered=True)

        modelflux, modelflux_ivar = sdsspy.dered_fluxes(st['extinction'], 
                                                        st['modelflux'], 
                                                        st['modelflux_ivar'])
        modelmag_dered, modelmag_dered_err = sdsspy.nmgy2mag(modelflux, modelflux_ivar)

        for f in sdsspy.FILTERCHARS:
            fnum = sdsspy.FILTERNUM[f]
            if 'devflux' in st.dtype.names:
                out_dict['cmodelflux_'+f] = cmodelmag_dered[:,fnum]
                out_dict['cmodelflux_ivar_'+f] = cmodelmag_dered_err[:,fnum]
                out_dict['cmodelmag_dered_'+f] = cmodelmag_dered[:,fnum]
                out_dict['cmodelmag_dered_err_'+f] = cmodelmag_dered_err[:,fnum]
            out_dict['modelmag_dered_'+f] = modelmag_dered[:,fnum]
            out_dict['modelmag_dered_err_'+f] = modelmag_dered_err[:,fnum]

        out_dict['photoid'] = sdsspy.photoid(st)

        self.set_maskflags(out_dict)

        # we will add an "survey_primary" column if we are not explicitly 
        # selecting survey_primary
        if self.type not in ['primgal','primstar']:
            survey_primary = numpy.zeros(st.size, dtype='i1')
            w=self.get_primary_indices(st)
            if w.size > 0:
                survey_primary[w] = 1
            out_dict['survey_primary'] = survey_primary

        return out_dict


    def get_primary_indices(self, objs, run_primary=False):
        if run_primary:
            primary = self.flags.val('resolve_status','run_primary')
        else:
            primary = self.flags.val('resolve_status','survey_primary')

        w=where1((objs['resolve_status'] & primary) != 0)
        return w

    def add_cmodelflux(self):
        """
        I forgot to add errors...
        """
        c = self.open_columns()
        for filt in ['u','g','r','i','z']:
            print("filter:",filt)
            fdev         = c['fracpsf_'+filt][:]
            devflux      = c['devflux_'+filt][:]
            devflux_ivar = c['devflux_ivar_'+filt][:]
            expflux      = c['expflux_'+filt][:]
            expflux_ivar = c['expflux_ivar_'+filt][:]
            flux, ivar   = sdsspy.util._make_cmodelflux_1band(fdev,
                                                              devflux, devflux_ivar,
                                                              expflux, expflux_ivar)
            c.write_column('cmodelflux_'+filt, flux, create=True)
            c.write_column('cmodelflux_ivar_'+filt, ivar, create=True)


    def add_dered_err_to_columns(self):
        """
        I forgot to add errors...
        """
        c = self.open_columns()

        for filt in ['u','g','r','i','z']:
            print("reading extinction for",filt)
            ext = c['extinction_'+filt][:]
            for type in ['cmodelmag','modelmag']:
                outtag = '%s_dered_err_%s' % (type,filt)
                print("making:",outtag)

                print("  reading fluxes")
                if type == 'cmodelmag':
                    fdev         = c['fracpsf_'+filt][:]
                    devflux      = c['devflux_'+filt][:]
                    devflux_ivar = c['devflux_ivar_'+filt][:]
                    expflux      = c['expflux_'+filt][:]
                    expflux_ivar = c['expflux_ivar_'+filt][:]
                    print("    making cmodelmag")
                    flux, ivar   = sdsspy.util._make_cmodelflux_1band(fdev,
                                                                      devflux, devflux_ivar,
                                                                      expflux, expflux_ivar)
                else:
                    flux = c['modelflux_'+filt][:]
                    ivar = c['modelflux_ivar_'+filt][:]

                print("  dereddening")
                flux, ivar = sdsspy.dered_fluxes(ext, flux, ivar)
                print("  making mags")
                mag, magerr = sdsspy.nmgy2mag(flux, ivar)

                print("  writing data to",outtag)
                c.write_column(outtag, magerr, create=True)

            


    def create_indices(self):
        c = self.open_columns()

        cnames = ['cmodelmag_dered_r','modelmag_dered_r',
                  'run','camcol','ra','dec','thing_id','photoid',
                  'inbasic','ingood','intycho']
        if self.type not in ['primgal','primstar']:
            cnames += ['survey_primary']

        for n in cnames:
            if n in c:
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
        #dirname += '_new'
        dir=os.environ['SWEEP_REDUCE']
        dir = path_join(dir,dirname)
        if not os.path.exists(dir):
            os.makedirs(dir)
                        
        dir=path_join(dir ,self.type+'.cols')
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




