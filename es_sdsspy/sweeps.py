import os,sys
from sys import stdout

import numpy
from numpy import where
import sdsspy
import esutil
from esutil.ostools import getenv_check, path_join

import datetime
import pprint

import columns

flagdict = {'read_failed':2**0,
            'noresult': 2**1,
            'nofile': 2**2,
            'proc_failed': 2**3}



def open_columns(type, primary=True):
    sel = ColumnSelector(type, primary=primary)
    dir = sel.columns_dir()
    return columns.Columns(dir)

def verify(proctype,
           procrun, 
           require_result=False,
           runs=None, 
           camcols=None, 
           reload=False):

    p=Proc(proctype,procrun,'junk')
    flag_status = p.verify(runs=runs,camcols=camcols,reload=reload)
    return flag_status

class Proc():
    def __init__(self, 
                 proctype, 
                 procrun, 
                 type='gal', 
                 nper=1, 
                 minscore=0.1,
                 reload=False, 
                 allow_crash=False):

        self.allow_crash=allow_crash
        self.init(proctype,procrun,type,nper,minscore,reload=reload)

    def init(self, 
             proctype, 
             procrun, 
             type='gal', 
             nper=1, 
             minscore=0.1,
             reload=False):
        self.proctype=proctype
        self.procrun=procrun
        self.type=type
        self.minscore=minscore


        self.nper = nper


    def load_runlist(self, reload=False):
        if not hasattr(Proc, '_runs') or reload:
            stdout.write("Getting run list, minscore=%s\n" % self.minscore)
            win = sdsspy.window.Window()
            Proc._runs, Proc._reruns = win.runlist(self.minscore)

    
    def process_runs(self, pyclass, runs=None, camcols=None, **keys):
        self.load_runlist()
        self.split_runlist()

        resdir = getenv_check('PHOTO_RESOLVE')
        calibdir = getenv_check('PHOTO_CALIB')
        sweepdir = getenv_check('PHOTO_SWEEP')
        redux = getenv_check('PHOTO_REDUX')

        runs2use,reruns2use = self.matchruns(runs)
        if camcols is None:
            camcols=[1,2,3,4,5,6]
        camcols2use = numpy.array(camcols, ndmin=1, copy=False) 

        for i in range(len(runs2use)):
            run = runs2use[i]
            rerun=reruns2use[i]

            for camcol in camcols2use:
                stdout.write("\nProcessing run: %s rerun: %s "
                             "camcol: %s\n\n" % (run,rerun,camcol))
                output_file = \
                    self.output_file(run=run,rerun=rerun,camcol=camcol)
                #stdout.write("Will output to file: %s\n" % output_file)
                # begin the status info
                now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
                status = \
                    {'run':run,
                     'rerun':rerun,
                     'photo_resolve':resdir,
                     'photo_calib':calibdir,
                     'photo_sweep':sweepdir,
                     'photo_redux':redux,
                     'esutilvers': esutil.version(),
                     'date': now,
                     'output_file':output_file}

                res,flags = self.process_sweep(run,rerun,camcol,pyclass,**keys)
                status['flags'] = flags
                if res is None:
                    nres = 0
                else:
                    nres = res.size
                status['nresult'] = nres

                print 'status:'
                pprint.pprint(status)
                self.write_result(output_file, status, res)
                del res


    def write_result(self, output_file, status, result):
        dir=os.path.dirname(output_file)
        if not os.path.exists(dir):
            stdout.write("Creating output directory: %s\n" % dir)
            os.makedirs(dir)

        if result is None:
            result = numpy.array([])

        esutil.io.write(output_file, result, 
                        header=status,
                        verbose=True,
                        clobber=True)

    def process_sweep(self,run,rerun,camcol,pyclass,**keys):
        flags = 0
        
        # process() method must take the objs array and **keys. 
        # on success, return a numpy.ndarray or superclass thereof.

        try:
            objs=self.read_sweep(run,rerun,camcol)
        except:
            info=sys.exc_info()
            stdout.write('Error: %s\n' % info[1])
            flags |= flagdict['read_failed']
            return None,flags

        # construct a new object
        pyobj = pyclass(objs, **keys)
        if self.allow_crash:
            result = pyobj.process(**keys)
        else:
            try:
                result = pyobj.process(**keys)
            except:
                info=sys.exc_info()
                stdout.write('Error: %s\n' % info[1])
                flags |= flagdict['proc_failed']
                return None,flags

        if result is None:
            flags |= flagdict['noresult']

        del objs
        return result, flags

    def columns_dir(self):
        dir=self.output_dir()
        coldir = dir+'.cols'
        return coldir

    def output_dir(self):
        dir=getenv_check('SWEEP_REDUCE')
        dir = path_join(dir,self.proctype,self.procrun)
        return dir

    def output_file(self, 
                    run=None, 
                    rerun=None, 
                    camcol=None, 
                    gather=False, 
                    collate=False,
                    extra=None):

        dir = self.output_dir()
        fname=[self.proctype,self.type,self.procrun]

        if not gather and not collate:
            if run is None or rerun is None or camcol is None:
                raise ValueError("if not gather and not collate, send "
                                 "run,rerun,camcol")
            fname.append(sdsspy.files.run2string(run))
            fname.append(sdsspy.files.camcol2string(camcol))
            fname.append(sdsspy.files.rerun2string(rerun))
        elif collate:
            fname.append('collate')
        elif gather:
            fname.append('gather')

        if extra is not None:
            fname.append(extra)

        #fname = '-'.join(fname) +'.fits'
        fname = '-'.join(fname) +'.rec'
        fname = os.path.join(dir, fname)
        return fname

    def read_sweep(self, run, rerun, camcol, type=None):
        if type is None:
            type=self.type
        data = sdsspy.files.read('calibObj.'+type,run,camcol,rerun=rerun,
                                 verbose=True,lower=True)
        return data

                         
    def split_runlist(self):
        self.split_runs = esutil.numpy_util.splitarray(self.nper, Proc._runs)
        self.split_reruns = esutil.numpy_util.splitarray(self.nper, Proc._reruns)

    def matchruns(self,runs=None, reload=False):
        self.load_runlist(reload=reload)

        if runs is None:
            return Proc._runs, Proc._reruns

        uruns = numpy.unique1d(runs)

        m1, m2 = esutil.numpy_util.match(uruns, Proc._runs)

        if m1[0] != -1:
            match_runs = Proc._runs[m2]
            match_reruns = Proc._reruns[m2]
            return match_runs, match_reruns
        else:
            # just return arrays [-1]
            return m1,m2


    def verify(self, 
               runs=None, 
               camcols=None, 
               require_result=False,
               reload=False):

        if camcols is None:
            camcols=[1,2,3,4,5,6]
        runs2use,reruns2use = self.matchruns(runs, reload=reload)
        camcols2use = numpy.array(camcols, ndmin=1, copy=False) 

        stdout.write("Verifying %s runs\n" % len(runs2use))


        nrun = len(runs2use)
        ncamcol=len(camcols)

        ntot = nrun*ncamcol

        stdout.write("%3d/%3d" % (0,nrun))

        flags = numpy.zeros(nrun*6, 'i4')

        ii=0
        for irun in range(nrun):
            if (irun+1) % 10 == 0:
                stdout.write('\b'*7)
                stdout.write("%3d/%3d" % ((irun+1),nrun) )
                stdout.flush()

            run=runs2use[irun]
            rerun=reruns2use[irun]

            for camcol in camcols:
                output_file = \
                    self.output_file(run=run,rerun=rerun,camcol=camcol)

                if not os.path.exists(output_file):
                    err="missing file: %s.\n    Halting" % output_file
                    raise RuntimeError(err)
                
                status=esutil.sfile.read_header(output_file)
                flags[ii] = status['flags']
                ii += 1

        stdout.write('\b'*7)
        stdout.write("%3d/%3d" % (nrun,nrun))
        stdout.write('\n')

        stdout.write("Checking flags\n")
        flag_status={}
        for flagname in flagdict:
            flag=flagdict[flagname]
            wbad, = where( (flags & flag) != 0)
            if wbad.size > 0:
                stdout.write("    %s: %3d/%3d\n" % (flagname,wbad.size,ntot))

            flag_status[flagname] = wbad.size
        
        for flagname in flag_status:
            if flag_status[flagname] != 0:
                if flagname == 'noresult' and not require_result:
                    stdout.write("No results required, ignoring noresult flag\n")
                else:
                    raise ValueError("Bad files found for flag : '%s'" % flagname)
        return flag_status


    def make_columns(self):
        """
        This is for the specialized runs.

        If you want everything from all sweeps or just primary, use ColumnSelector.
        """
        import columns

        coldir = self.columns_dir()
        if os.path.exists(coldir):
            raise ValueError("coldir already exists: %s\n"
                             "\tPlease start fresh" % coldir)

        self.verify()

        c=columns.Columns(coldir)
        c.create()

        rerun=self._reruns[0]
        irun=1
        nrun=len(self._runs)
        for run in sorted(self._runs):
            stdout.write("%s/%s\n" % (irun,nrun))
            irun += 1

            for camcol in [1,2,3,4,5,6]:
                f=self.output_file(run=run,rerun=rerun,camcol=camcol)

                status = esutil.sfile.read_header(f)

                if status['flags'] == 0:
                    data = esutil.io.read(f)
                    c.write_columns(data)
                else:
                    stdout.write("no results in file %s\n" % f)

                del status


class ColumnSelector:
    """

    Put all the objects into a columns database.  It will go under 
        ~/sweeps_reduce/{PHOTO_SWEEP}/{type}.cols
    If doing primary, it will be
        ~/sweeps_reduce/{PHOTO_SWEEP}/prim{type}.cols

    More complex selections can easily be added.

    """
    def __init__(self, type, minscore=0.1, primary=True):
        self.type=type
        self.minscore=minscore
        self.primary=primary

        self.flags = sdsspy.flags.Flags()
        self.masktypes = ['basic','good','tycho']

    def make_indices(self, columns=None, tempdir='/dev/shm'):
        if columns is None:
            columns=['run','camcol','field','id','photoid','thing_id',
                     'ra','dec','htmid07','htmid10','inbasic','ingood','intycho']
        self.open_columns()

        for c in columns:
            self.cols[c].create_index(force=True, verbose=True, tempdir=tempdir)

    def make_meta(self):
        import datetime
        date = str( datetime.datetime.now() )
        psweep_path = os.getenv('PHOTO_SWEEP')
        psweep = os.path.basename(psweep_path)
        out = {}
        out['date'] = date
        out['photo_sweep'] = psweep
        out['photo_sweep_path'] = psweep_path
        out['primary'] = self.primary
        out['minscore'] = self.minscore
        out['type'] = self.type

        return out

    def write_meta(self):
        meta = self.make_meta()
        self.cols.write_column('meta', meta)


    def process_all(self):
        self.open_columns_new()
        self.write_meta()
        return

        win = sdsspy.window.Window()
        runs, reruns = win.runlist(self.minscore)
        nrun=len(runs)

        for i in xrange(nrun):
            stdout.write("\nRun: %s/%s\n" % ((i+1),nrun))
            self.process_run(runs[i], reruns[i])

    def process_run(self, run, rerun):
        for camcol in [1,2,3,4,5,6]:
            data=sdsspy.files.read('calibobj',run,camcol=camcol,
                                   rerun=rerun, type=self.type,
                                   verbose=True)
            output = self.select(data)
            if output is not None:
                self.cols.write_columns(output)

    def open_columns_new(self):
        coldir=self.columns_dir()
        if os.path.exists(coldir):
            raise ValueError("coldir exists: '%s'. Start fresh" % coldir)

        stdout.write("Opening columns: %s\n" % coldir)
        self.cols = columns.Columns(coldir)
        self.cols.create()

    def open_columns(self):
        coldir=self.columns_dir()
        stdout.write("Opening columns: %s\n" % coldir)
        self.cols = columns.Columns(coldir)


    def name(self):
        if self.primary:
            name='prim'
        else:
            name=''
        return '{name}{type}'.format(name=name,type=self.type)

    def columns_dir(self):
        if 'SWEEP_REDUCE' not in os.environ:
            raise ValueError("SWEEP_REDUCE is not set")
        dir = os.environ['SWEEP_REDUCE']

        if 'PHOTO_SWEEP' not in os.environ:
            raise ValueError("PHOTO_SWEEP is not set")
        sweep=os.environ['PHOTO_SWEEP'] 
        sweep = os.path.basename(sweep)

        dir = path_join(dir, sweep, self.name()+'.cols')
        return dir

    def keep_dtype(self, orig_dtype):
        dtype = list(orig_dtype)
        extra=[('photoid','i8'),
               ('modelmag_dered', 'f4', 5),
               ('modelmag_dered_err', 'f4', 5),
               ('cmodelmag_dered', 'f4', 5),
               ('cmodelmag_dered_err', 'f4', 5),
               ('htmid07','i4'),
               ('htmid10','i4'),
               ('inbasic','i1'),
               ('ingood','i1'),
               ('intycho','i1')]
        dtype += extra
        return dtype

        dtype = \
            [('photoid', 'i8'),
             ('run', 'i2'),
             ('rerun', '|S3'),
             ('camcol', 'i1'),
             ('field', 'i2'),
             ('id', 'i2'),
             ('thing_id','i4'),
             ('ra', 'f8'),
             ('dec', 'f8'),
             ('flags', 'i4', 5),
             ('flags2', 'i4', 5),
             ('objc_flags', 'i4'),
             ('objc_flags2', 'i4'),
             ('psfflux', 'f4', 5),
             ('psfflux_ivar', 'f4', 5),
             ('modelflux', 'f4', 5),
             ('modelflux_ivar', 'f4', 5),
             ('modelmag_dered', 'f4', 5),
             ('modelmag_dered_err', 'f4', 5),
             ('extinction', 'f4', 5),
             ('psf_fwhm', 'f4', 5),
             ('calib_status','i4',5),
             ('htmid07','i4'),
             ('htmid10','i4'),
             ('inbasic','i1'),
             ('ingood','i1'),
             ('intycho','i1')]
        return dtype

    def select(self, objs): 
        """
        Select the primary objects and create the output structure
        """
        if self.primary:
            stdout.write("Getting primary\n")
            logic = self.resolve_logic(objs)
            w, = numpy.where(logic)

            stdout.write('Keeping %s/%s\n' % (w.size,objs.size))
            if w.size > 0:
                result = self.extract_fields(objs[w])
            else:
                result = None
        else:
            stdout.write('nobj: %s\n' % objs.size)
            result = self.extract_fields(objs)

        return result


    def resolve_logic(self, objs):
        primary_flag = self.flags.val('resolve_status','survey_primary')
        return (objs['resolve_status'] & primary_flag) != 0

    def extract_fields(self, objs):
        """
        Extract tags from the input structure and add some if needed.
        """

        dtype=self.keep_dtype(objs.dtype.descr)
        new = numpy.zeros(objs.size, dtype=dtype)
        esutil.numpy_util.copy_fields(objs, new)
    
        # calculate new fields
        new['photoid'] = sdsspy.photoid(new)

        flux,ivar=new['modelflux'],new['modelflux_ivar']
        ext=new['extinction']
        flux_dered,ivar_dered = sdsspy.dered_fluxes(ext, flux, ivar=ivar)


        mag_dered,err_dered = sdsspy.util.nmgy2mag(flux_dered, ivar=ivar_dered)
        new['modelmag_dered'] = mag_dered
        new['modelmag_dered_err'] = err_dered

        cmag_dered,cmagerr_dered = sdsspy.make_cmodelmag(objs, dered=True)
        new['cmodelmag_dered'] = cmag_dered
        new['cmodelmag_dered_err'] = cmagerr_dered
        
        self.set_htmids(new)
        self.set_maskflags(new)

        return new

    def set_htmids(self, new):
        # add htmid
        ra = new['ra']
        dec = new['dec']
        stdout.write("    Getting htm depth ")
        for depth in [7,10]:
            stdout.write('%s ' % depth)
            name='htmid%02i' % depth
            h=esutil.htm.HTM(depth)
            new[name] = h.lookup_id(ra, dec)
        stdout.write('\n')
    def set_maskflags(self, new):
        self.load_masks()
        ra = new['ra']
        dec = new['dec']
        stdout.write("    Checking against mask: ")
        for masktype in self.masktypes:
            name = 'in'+masktype
            stdout.write('%s ' % name)

            maskflags = numpy.zeros(ra.size, dtype='i1')


            new[name] = self.masks[masktype].Contains(ra,dec,'eq')
     
        stdout.write('\n')

        return new


    def load_masks(self):
        if not hasattr(self,'masks'):
            self.masks = {}
            for masktype in self.masktypes:
                self.masks[masktype] = sdsspy.stomp_maps.load('boss',masktype, 
                                                              verbose=True)




