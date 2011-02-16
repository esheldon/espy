"""
First create the pbs files, e.g. with a function like
    pbs_gal04
Submit using the .sh file

That will also create the python script to load the columns.
    loadcol-sdssgal-prim04-byrun.py

Then add the spatial information: HTM and stomp maps
    add_spatial_info(procrun)

"""
import esutil
from esutil.ostools import getenv_check, path_join, expand_path
import os
import sdsspy
import numpy
from numpy import where
import sys
from sys import stdout
import time

import pprint

def pbs_gal04(primonly=True):
    proctype='sdssgal'
    if primonly:
        procrun='prim04'
    else:
        procrun='04'
    gs = GalSelector(primonly=primonly)
    gs.create_pbs_byrun(proctype,procrun)

def pbs_gal03(primonly=False):
    proctype='sdssgal'
    if primonly:
        procrun='prim03'
    else:
        procrun='03'
    gs = GalSelector(primonly=primonly)
    gs.create_pbs_byrun(proctype,procrun)

class GalSelector():
    """
    Note for primonly, the only cuts made are on resolve.
    No mag, flags, or other cuts are made.
    """
    def __init__(self, objs=None, **keys):
        self.init(objs, **keys)

    def init(self, objs=None, **keys):
        self.objs = objs
        self.max_rmag = keys.get('max_rmag', 22.0)
        self.primonly=keys.get('primonly',False)

        self.flags = sdsspy.flags.Flags()

    def process(self, **keys):
        stdout.write("Getting logic\n")
        logic = self.select_logic()
        stdout.write("Selecting good galaxies\n")
        w, = numpy.where(logic)
        if w.size > 0:
            result = self.objs[w]

            result = self.extract_fields(result)
        else:
            result = None
        stdout.write('Keeping %s/%s\n' % (w.size,self.objs.size))
        return result

    def keep_dtype(self):
        dtype = \
            [('photoid', 'i8'),
             ('thing_id','i4'),  # will be -1 for some non-primary objects
             ('run', 'i2'),
             ('rerun', '|S3'),
             ('camcol', 'i1'),
             ('field', 'i2'),
             ('id', 'i2'),
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
             ('cmodelmag_dered', 'f4', 5),
             ('cmodelmag_dered_err', 'f4', 5),
             ('extinction', 'f4', 5),
             ('psf_fwhm', '>f4', 5)]

        if self.primonly:
            dtype.append( ('calib_status','i4',5) )

        return dtype

    def extract_fields(self, objs):
        """
        Extract tags from the input structure and add some if needed.
        """

        dtype=self.keep_dtype()
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

        return new


    def select_logic(self):
        if self.primonly:
            logic = self.resolve_logic()
        else:
            mag_logic = self.mag_logic()
            resolve_logic = self.resolve_logic()
            flag_logic = self.flag_logic()

            logic = mag_logic & resolve_logic & flag_logic

        return logic



    def mag_logic(self):
        cmodel_dered = sdsspy.make_cmodelmag(self.objs, doerr=False, dered=True)
        mag_logic = cmodel_dered[:,2] < self.max_rmag
        return mag_logic

    def resolve_logic(self):
        primary_flag = self.flags.val('resolve_status','survey_primary')
        return (self.objs['resolve_status'] & primary_flag) != 0





        

    def flag_logic(self):

        # make sure detected in r and i in binned1,2, or 4
        binned_logic = self.binned_logic()

        # make sure no non-photometric data
        calib_logic = self.calib_logic()

        # checking lots of object1 flags
        object1_logic = self.object1_logic()

        # combined logic
        logic = binned_logic & calib_logic & object1_logic

        return logic

    def object1_logic(self):
        satur = self.flags.val('object1','satur')
        bright = self.flags.val('object1','bright')
        too_many_peaks = self.flags.val('object1','deblend_too_many_peaks')

        blended = self.flags.val('object1','blended')
        nodeblend = self.flags.val('object1','nodeblend')

        peakcenter = self.flags.val('object1','peakcenter')
        notchecked = self.flags.val('object1','notchecked')
        noprofile = self.flags.val('object1','noprofile')

        oflags = self.objs['objc_flags']

        # famously worded as double negatives
        bdb_logic = \
            ( (oflags & blended) == 0) | ((oflags & nodeblend) != 0)

        # now combine logic
        logic = \
            ((oflags & satur) == 0) \
            & ((oflags & bright) == 0) \
            & ((oflags & too_many_peaks) == 0) \
            & ((oflags & peakcenter) == 0) \
            & ((oflags & notchecked) == 0) \
            & ((oflags & noprofile) == 0) \
            & bdb_logic

        return logic


    def binned_logic(self):
        binned1 = self.flags.val('object1','binned1')
        binned2 = self.flags.val('object1','binned2')
        binned4 = self.flags.val('object1','binned4')

        rflags = self.objs['flags'][:,2]
        iflags = self.objs['flags'][:,3]

        r_binned_logic = \
            ((rflags & binned1) != 0 )  \
            | ((rflags & binned2) != 0 )  \
            | ((rflags & binned4) != 0 )

        i_binned_logic = \
            ((iflags & binned1) != 0 )  \
            | ((iflags & binned2) != 0 )  \
            | ((iflags & binned4) != 0 )

        return r_binned_logic & i_binned_logic




    def calib_logic(self, objs=None):
        if objs is None:
            objs = self.objs
        calib_status = objs['calib_status']
        flagval = self.flags.val('calib_status', 'photometric')
        logic = \
            ((calib_status[:,0] & flagval) != 0)   \
            & ((calib_status[:,1] & flagval) != 0) \
            & ((calib_status[:,2] & flagval) != 0) \
            & ((calib_status[:,3] & flagval) != 0) \
            & ((calib_status[:,4] & flagval) != 0)


        return logic
        

    def create_pbs_byrun(self,
                         proctype,
                         procrun,
                         espy_v=None,
                         esutil_v=None,
                         stomp_v=None,
                         photo_sweep=None,
                         photo_calib=None,
                         photo_resolve=None):

        import sdsspy
        import pbs

        p = sdsspy.sweeps.Proc(proctype,procrun)

        extra=''
        if self.primonly:
            extra=', primonly=True'

        pbsdir=path_join('~/pbs',proctype,procrun)
        pbsdir = expand_path(pbsdir)

        if not os.path.exists(pbsdir):
            stdout.write("Making directory: %s\n" % pbsdir)
            os.makedirs(pbsdir)

        if espy_v is None: 
            espy_v = "-r ~esheldon/exports/espy-work"
        if esutil_v is None:
            esutil_v="-r ~esheldon/exports/esutil-work"
        if stomp_v is None:
            stomp_v="-r ~esheldon/exports/stomp-work"

        if photo_sweep is None:
            photo_sweep=os.getenv("PHOTO_SWEEP")
        if photo_resolve is None:
            photo_resolve=os.getenv("PHOTO_RESOLVE")
        if photo_calib is None:
            photo_calib=os.getenv("PHOTO_CALIB")

        setups = \
            [espy_v,esutil_v,stomp_v]
        setups = ['setup '+s for s in setups]

        setups += ['export PHOTO_SWEEP=%s' % photo_sweep,
                   'export PHOTO_CALIB=%s' % photo_calib,
                   'export PHOTO_RESOLVE=%s' % photo_resolve]
        
        
        qsub_filename=path_join(pbsdir, 'submit-%s-%s-byrun.sh' % (proctype,procrun))
        loadcol_filename=path_join(pbsdir, 'loadcol-%s-%s-byrun.py' % (proctype,procrun))

        loadcol_file=open(loadcol_filename,'w')
        loadcol_file.write("import sdsspy\n")
        loadcol_file.write("p=sdsspy.sweeps.Proc('%s','%s')\n" % (proctype,procrun))
        loadcol_file.write("p.make_columns()\n")
        loadcol_file.close()


        qsub_file=open(qsub_filename, 'w')

        irun=0
        nrun=len(p._runs)
        for run in sorted(p._runs):
            stdout.write("Creating pbs file for run: %s\n" % run)

            # create the pbs file name
            rstr=sdsspy.files.run2string(run)
            job_name = '%s-%s' % (proctype,rstr)

            pbsfilename='%s-%s-%s.pbs' % (proctype, procrun, rstr)
            pbsfilename=path_join(pbsdir, pbsfilename)

            # the python commands to execute
            python_commands = """
import sdsspy
import sdssgal
gs=sdssgal.GalSelector
p=sdsspy.sweeps.Proc('sdssgal','%s','gal')
p.process_runs(gs, runs=%s%s)""" % (procrun,run,extra)

            ppy = pbs.PBSPython(pbsfilename,python_commands,
                                job_name=job_name, setups=setups)
            ppy.write()

 
            # write to the submission file
            qsub_file.write('echo -n "%s/%s " \n' % (irun+1,nrun))
            qsub_file.write('qsub %s\n' % pbsfilename)


            irun += 1

        qsub_file.close()


def calib_logic(calib_status):

    flagval = sdsspy.flagval('calib_status', 'photometric')
    logic = \
        ((calib_status[:,0] & flagval) != 0)   \
        & ((calib_status[:,1] & flagval) != 0) \
        & ((calib_status[:,2] & flagval) != 0) \
        & ((calib_status[:,3] & flagval) != 0) \
        & ((calib_status[:,4] & flagval) != 0)

    return logic


def resolve_logic(resolve_status):
    primary_flag = sdsspy.flagval('resolve_status','survey_primary')
    return (resolve_status & primary_flag) != 0


def binned_logic_ri(flags):
    binned1 = sdsspy.flagval('object1','binned1')
    binned2 = sdsspy.flagval('object1','binned2')
    binned4 = sdsspy.flagval('object1','binned4')

    rflags = flags[:,2]
    iflags = flags[:,3]

    r_binned_logic = \
        ((rflags & binned1) != 0 )  \
        | ((rflags & binned2) != 0 )  \
        | ((rflags & binned4) != 0 )

    i_binned_logic = \
        ((iflags & binned1) != 0 )  \
        | ((iflags & binned2) != 0 )  \
        | ((iflags & binned4) != 0 )

    return r_binned_logic & i_binned_logic

def object1_logic(objc_flags):
    satur = sdsspy.flagval('object1','satur')
    bright = sdsspy.flagval('object1','bright')
    too_many_peaks = sdsspy.flagval('object1','deblend_too_many_peaks')

    blended = sdsspy.flagval('object1','blended')
    nodeblend = sdsspy.flagval('object1','nodeblend')

    peakcenter = sdsspy.flagval('object1','peakcenter')
    notchecked = sdsspy.flagval('object1','notchecked')
    noprofile = sdsspy.flagval('object1','noprofile')


    # famously worded as double negatives
    bdb_logic = \
        ( (objc_flags & blended) == 0) | ((objc_flags & nodeblend) != 0)

    # now combine logic
    logic = \
        ((objc_flags & satur) == 0) \
        & ((objc_flags & bright) == 0) \
        & ((objc_flags & too_many_peaks) == 0) \
        & ((objc_flags & peakcenter) == 0) \
        & ((objc_flags & notchecked) == 0) \
        & ((objc_flags & noprofile) == 0) \
        & bdb_logic

    return logic






# this is  just for testing
def run_sweeproc(runs=None, camcols=[1,2,3,4,5,6], allow_crash=False):
    gs = GalSelector

    p = sdsspy.sweeps.Proc('sdssgal','gal01','gal', allow_crash=allow_crash)

    p.process_runs(gs, runs=runs, camcols=camcols)




def open_columns(procrun, verbose=False):
    import columns
    coldir=columns_dir(procrun)
    stdout.write("Opening galaxies columns dir: '%s'\n" % coldir)
    if not os.path.exists(coldir):
        raise ValueError("missing coldir: %s" % coldir)
    return columns.Columns(coldir, verbose=verbose)

def columns_dir(procrun):
    sweep_reduce=getenv_check('SWEEP_REDUCE')
    coldir=path_join(sweep_reduce,'sdssgal','%s.cols' % procrun)
    return coldir

def add_spatial_info(procrun):
    import stomp

    c=open_columns(procrun, verbose=True)


    stdout.write("Reading ra/dec\n")
    ra=c['ra'][:]
    dec=c['dec'][:]


    for depth in [7,10]:
        tm0=time.time() 
        stdout.write("Getting HTM id at depth %d\n" % depth)
        h=esutil.htm.HTM(depth)
        htmid = h.lookup_id(ra, dec)
        tm1=time.time()
        esutil.misc.ptime(tm1-tm0)

        colname='htmid%02i' % depth
        stdout.write("Creating htm id column '%s'\n" % colname)
        c.write_column(colname, htmid, meta={'depth':depth},create=True)


    stdout.write("Getting window info\n")

    for masktype in ['basic','good','tycho']:
        colname = 'in'+masktype

        maskflags = numpy.zeros(ra.size, dtype='i1')

        m = sdsspy.stomp_maps.load('boss',masktype)



        stdout.write("Checking against mask\n"); stdout.flush()
        tm0=time.time() 
        maskflags[:] = m.Contains(ra,dec,'eq')
        tm1=time.time() 
        esutil.misc.ptime(tm1-tm0)
        stdout.flush()

        c.write_column(colname, maskflags, create=True)

        del m
        del maskflags
                    

def plot_masked(procrun, show=False):
    import biggles

    c=open_columns(procrun)

    sz=c['inbasic'].size
    nrand = sz/100
    nrand = sz/1000
    stdout.write('getting random subset %s/%s\n' % (nrand,sz))

    r=esutil.numpy_util.random_subset(sz, nrand)

    #inbasic = c['inbasic'][r]
    #wbasic,=where(inbasic == 1)
    #del inbasic
    intycho=c['intycho'][r]
    wtycho,=where(intycho == 1)
    del intycho
    
    ra=c['ra'][r]
    ra = esutil.coords.shiftlon(ra, 90.0)
    #ra = esutil.coords.shiftlon(ra)
    dec=c['dec'][r]

    plt=biggles.FramedPlot()

    pall = biggles.Points( ra, dec, type='dot' )
    ptycho = biggles.Points( ra[wtycho], dec[wtycho], type='dot', color='red' )

    # for the legend
    lall = biggles.Points( -1000, -1000, type='filled circle')
    lall.label = 'All'
    ltycho = biggles.Points( -1000, -1000, type='filled circle', color='red')
    ltycho.label = 'Within Tycho'

    plt.xlabel = 'RA-90'
    plt.ylabel = 'DEC'

    plt.add(pall)
    plt.add(ptycho)
    plt.add( biggles.PlotKey(0.6,0.9, [lall,ltycho] ) )

    if show:
        plt.show()
    #png='/home/users/esheldon/www/tmp/plots/sdssgal-%s-boss-mask.png' % procrun
    #stdout.write("writing plot to '%s'\n" % png)
    #plt.write_img(800,800,png)
    plotf='/home/users/esheldon/www/tmp/plots/sdssgal-%s-boss-mask.eps' % procrun
    stdout.write("writing plot to '%s'\n" % plotf)
    plt.write_eps(plotf)
    dpi=150
    esutil.ostools.exec_process('converter -d %s %s' % (dpi,plotf),verbose=True)

def plot_masked_progressive(procrun, show=False, numcut=100):
    import biggles

    c=open_columns(procrun)

    sz=c['inbasic'].size
    nrand = sz/numcut
    #nrand = sz/1000
    stdout.write('getting random subset %s/%s\n' % (nrand,sz))

    r=esutil.numpy_util.random_subset(sz, nrand)

    #inbasic = c['inbasic'][r]
    #wbasic,=where(inbasic == 1)
    #del inbasic
    intycho=c['intycho'][r]
    wtycho,=where(intycho == 1)

    calib_status = c['calib_status'][r]

    clogic = calib_logic(calib_status)
    wcalib, = where( (intycho == 1) & (clogic) )

    del intycho
    del clogic
    
    ra=c['ra'][r]
    ra = esutil.coords.shiftlon(ra, 90.0)
    #ra = esutil.coords.shiftlon(ra)
    dec=c['dec'][r]

    plt=biggles.FramedPlot()

    symsize=0.7
    pall = biggles.Points( ra, dec, type='dot',size=symsize )
    ptycho = biggles.Points( ra[wtycho], dec[wtycho], type='dot', color='blue',
                            size=symsize)
    pcalib = biggles.Points( ra[wcalib], dec[wcalib], type='dot', color='red',
                           size=symsize)

    # for the legend
    lall = biggles.Points( -1000, -1000, type='filled circle')
    lall.label = 'Primary'
    ltycho = biggles.Points( -1000, -1000, type='filled circle', color='blue')
    ltycho.label = 'Within Tycho'

    lcalib = biggles.Points( -1000, -1000, type='filled circle', color='red')
    lcalib.label = 'Tycho & photometric'

    plt.xlabel = 'RA-90'
    plt.ylabel = 'DEC'

    plt.add(pall)
    plt.add(ptycho)
    plt.add(pcalib)
    plt.add( biggles.PlotKey(0.6,0.9, [lall,ltycho,lcalib],size=1 ) )

    if show:
        plt.show()
    else:
        plotf='/home/users/esheldon/www/tmp/plots/sdssgal-%s-boss-mask-calib.eps' % procrun
        stdout.write("writing plot to '%s'\n" % plotf)
        plt.write_eps(plotf)
        dpi=150
        esutil.ostools.exec_process('converter -d %s %s' % (dpi,plotf),verbose=True)


   
    
def shiftra(ra, shift=None, wrap=False):
    """
    ranew = shiftra(ra, shift=90.0)
    ranew = shiftra(ra, wrap=True)
    """

    if wrap:
        ra_wrapped = ra.copy()
        w,=numpy.where( ra > 180 )
        if w.size > 0:
            ra_wrapped[w] -= 360
        
        return ra_wrapped

    if shift is not None:
        ra_shift = ra - shift
        w, = numpy.where(ra_shift < 0)
        if w.size > 0:
            ra_shift[w] += 360

        return ra_shift

    return ra






