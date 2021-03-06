"""

classes:

    RegaussSweep
    RegaussAtlas
    SweepCache

    Collator: 
        - Collate into columns

    Zphot:
        - add a match index into photoz column databases
        - add interpolated inverse critical density curves

    match_zphot: Match to the photoz outputs and add a column with this
        match index to the columns database.

    RunTester
    Tester

"""
from __future__ import print_function
import os, sys
import glob


import numpy
from numpy import where,sqrt
import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join, expand_path

import images

def open_columns(procrun, sweeptype):
    coll = Collator(procrun, sweeptype)
    return coll.open_columns()

def output_dir(procrun, sweeptype):
    coll = Collator(procrun, sweeptype)
    return coll.output_dir()


def zphot_match(procrun, pzrun, sweeptype):
    '''
    add a match index into photoz column databases

    The column name will be match_zphot{pzrun}
    '''
    import zphot

    print("opening regauss columns for procrun:",procrun)
    rgcols = open_columns(procrun, sweeptype)
    print("  #rows:",rgcols['photoid'].size)
    print("opening zphot columns for pzrun:",pzrun)
    pzcols = zphot.weighting.open_pofz_columns(pzrun)
    print("  #rows:",pzcols['photoid'].size)



    print("reading num from regauss")
    num = pzcols['num'][:]

    print("reading photoid from regauss")
    rg_photoid = rgcols['photoid'][:]
    print("reading photoid from zphot")
    pz_photoid = pzcols['photoid'][:]

    print("matching")
    mrg, mpz =  eu.numpy_util.match(rg_photoid, pz_photoid)
    print("  matched:",mrg.size)


    print("now determining which zphot are recoverable")
    w_recover = where1(num[mpz] > 0)
    print("  found:",w_recover.size)

    mrg = mrg[w_recover]
    mpz = mpz[w_recover]

    matches = numpy.empty(rg_photoid.size, dtype='i4')
    matches[:] = -1

    # must explicitly convert to i4
    matches[mrg] = numpy.array(mpz, dtype='i4')

    match_column = 'match_zphot%s' % pzrun
    print("Adding zphot match column:",match_column)
    rgcols.write_column(match_column, matches, create=True)
    print("Creating index")
    rgcols[match_column].create_index()

def make_scinv_chunk_scripts(procrun, scinv_sample, chunksize=1000000):
    print("using chunksize:",chunksize)

    outdir='~/lensing/proc/proc%s-add-scinv-%s' % (procrun,scinv_sample)
    outdir=os.path.expanduser(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cols = open_columns(procrun, 'gal')

    ntot = cols['run'].size
    nchunk = ntot/chunksize
    if (ntot % chunksize) > 0:
        nchunk += 1

    for chunk in xrange(nchunk):
        fname='chunk%03d.py' % chunk
        fname=os.path.join(outdir, fname)

        print(fname)
        with open(fname,'w') as fobj:
            prog="""import lensing
lensing.regauss.add_scinv_chunk('%s','%s',%s)\n"""
            
            prog=prog % (procrun, scinv_sample, chunk)

            fobj.write(prog)

def combine_scinv_chunks(procrun, scinv_sample, chunksize=1000000):
    from . import files, scat, sigmacrit
    import zphot

    print("using chunksize:",chunksize)

    scinv_conf = files.read_config('scinv', scinv_sample)
    cosmo = files.read_config('cosmo',scinv_conf['cosmo_sample'])
    pzrun=scinv_conf['pzrun']

    cols = open_columns(procrun, 'gal')

    ntot = cols['run'].size
    nchunk = ntot/chunksize
    if (ntot % chunksize) > 0:
        nchunk += 1

    corrfile=zphot.weighting.pofz_correction_file(pzrun)
    corrstruct = eu.io.read(corrfile)
    zs = (corrstruct['zmax']+corrstruct['zmin'])/2.
    corr = corrstruct['corr']

    scalc = sigmacrit.ScinvCalculator(scinv_conf['zlmin'],
                                      scinv_conf['zlmax'],
                                      scinv_conf['nzl'], 
                                      zs[0],
                                      zs[-1],
                                      H0=cosmo['H0'],
                                      omega_m=cosmo['omega_m'])


    zlvals = scalc.zlvals


    meta={'scinv_sample':scinv_sample,
          'zlmin':scinv_conf['zlmin'],
          'zlmax':scinv_conf['zlmax'],
          'nzl':scinv_conf['nzl'],
          'zlvals':zlvals}


    colname = scat.scinv_colname(scinv_sample)
    print("Writing to column:\n",colname)
    for chunk in xrange(nchunk):
        chunkname='%s-%03d' % (colname,chunk)

        print('  reading:',chunkname)
        data=cols[chunkname][:]

        print('  writing')
        if chunk == 0:
            cols.write_column(colname, data,meta=meta)
        else:
            cols.write_column(colname, data)



def add_scinv_chunk(procrun, scinv_sample, chunk, chunksize=1000000):
    """

    scinv are created using the *corrected* p(z).

    The nzl,zlmin,zlmax,zlvals will be in the meta data

    Procedure
    ---------
    - get the correction factor N(z)/sumpofz(z)
        calls zphot.weighting.pofz_correction(pzrun)

    - for each p(z) multiply by that factor
    - generate mean inverse critical density.


    """
    from . import files, scat, sigmacrit
    import zphot

    scinv_conf = files.read_config('scinv', scinv_sample)
    cosmo = files.read_config('cosmo',scinv_conf['cosmo_sample'])
    pzrun=scinv_conf['pzrun']


    cols = open_columns(procrun, 'gal')

    # the column which will hold the inverse critical density.
    # depending on keywords, we might want to raise an error
    colname = scat.scinv_colname(scinv_sample)
    colname='%s-%03d' % (colname,chunk)

    print("Writing to column:\n",colname)
    if colname in cols:
        if not clobber:
            raise ValueError("Column already exists")
        else:
            print("  removing column")
            cols[colname].delete()

    # get the matches
    zphot_matchname = scat.zphot_matchname(pzrun)

    if zphot_matchname not in cols:
        raise ValueError("zphot match column not "
                         "found: '%s'" % zphot_matchname)

    #
    # correction factor to apply to all p(z)
    #
    print("getting p(z) correction function\n")
    corrfile=zphot.weighting.pofz_correction_file(pzrun)
    corrstruct = eu.io.read(corrfile)
    zs = (corrstruct['zmax']+corrstruct['zmin'])/2.
    corr = corrstruct['corr']


    print("")
    scalc = sigmacrit.ScinvCalculator(scinv_conf['zlmin'],
                                      scinv_conf['zlmax'],
                                      scinv_conf['nzl'], 
                                      zs[0],
                                      zs[-1],
                                      H0=cosmo['H0'],
                                      omega_m=cosmo['omega_m'])

    zlvals = scalc.zlvals

    meta={'scinv_sample':scinv_sample,
          'zlmin':scinv_conf['zlmin'],
          'zlmax':scinv_conf['zlmax'],
          'nzl':scinv_conf['nzl'],
          'zlvals':zlvals}

    print("opening corresponding p(z) columns: '%s'\n" % pzrun)
    pzcols = zphot.weighting.open_pofz_columns(pzrun)

    # work on chunks
    print("using chunksize:",chunksize)
    ntot = cols['run'].size
    nchunk = ntot/chunksize
    if (ntot % chunksize) > 0:
        nchunk += 1

    print("nchunk:",nchunk)
    if chunk >= nchunk:
        raise ValueError("chunk should be in [0,%d]" % (nchunk-1,) )

    imin = chunk*chunksize
    imax = (chunk+1)*chunksize
    print("  %d %d (%d)" % (imin,imax,ntot))
    match = cols[zphot_matchname][imin:imax]
    
    nrows = match.size
    scinv = numpy.zeros((nrows, zlvals.size), dtype='f8') -9999.

    w,=numpy.where(match >= 0)
    print("    ",w.size,"with matches")
    if w.size > 0:
        # read the p(z) for the matches
        pofz = pzcols['pofz'][match[w]]

        pofz = pofz[:]*corr[:]

        for j in xrange(w.size):
            scinv[w[j], :] = scalc.calc_mean_scinv(pofz[j])
    # metadata only gets written once
    cols.write_column(colname, scinv, meta=meta)

def add_scinv(procrun, scinv_sample, 
              clobber=False,
              chunksize=100000):
    """

    Use the chunked versioin instead

    Add a column to the db that is the mean scinv as a function of lens
    redshift.

    scinv are created using the *corrected* p(z).

    The nzl,zlmin,zlmax,zlvals will be in the meta data

    Procedure
    ---------
    - get the correction factor N(z)/sumpofz(z)
        calls zphot.weighting.pofz_correction(pzrun)

    - for each p(z) multiply by that factor
    - generate mean inverse critical density.


    """
    from . import files, scat, sigmacrit
    import zphot

    scinv_conf = files.read_config('scinv', scinv_sample)
    cosmo = files.read_config('cosmo',scinv_conf['cosmo_sample'])
    pzrun=scinv_conf['pzrun']


    cols = open_columns(procrun, 'gal')

    # the column which will hold the inverse critical density.
    # depending on keywords, we might want to raise an error
    colname = scat.scinv_colname(scinv_sample)
    print("Writing to column:\n",colname)
    if colname in cols:
        if not clobber:
            raise ValueError("Column already exists")
        else:
            print("  removing column")
            cols[colname].delete()

    # get the matches
    zphot_matchname = scat.zphot_matchname(pzrun)

    if zphot_matchname not in cols:
        raise ValueError("zphot match column not "
                         "found: '%s'" % zphot_matchname)

    #
    # correction factor to apply to all p(z)
    #
    print("getting p(z) correction function\n")
    corrfile=zphot.weighting.pofz_correction_file(pzrun)
    corrstruct = eu.io.read(corrfile)
    zs = (corrstruct['zmax']+corrstruct['zmin'])/2.
    corr = corrstruct['corr']


    print("")
    scalc = sigmacrit.ScinvCalculator(scinv_conf['zlmin'],
                                      scinv_conf['zlmax'],
                                      scinv_conf['nzl'], 
                                      zs[0],
                                      zs[-1],
                                      H0=cosmo['H0'],
                                      omega_m=cosmo['omega_m'])


    zlvals = scalc.zlvals

    meta={'scinv_sample':scinv_sample,
          'zlmin':scinv_conf['zlmin'],
          'zlmax':scinv_conf['zlmax'],
          'nzl':scinv_conf['nzl'],
          'zlvals':zlvals}

    print("opening corresponding p(z) columns: '%s'\n" % pzrun)
    pzcols = zphot.weighting.open_pofz_columns(pzrun)

    # work on chunks
    print("using chunksize:",chunksize)
    ntot = cols['run'].size
    nchunk = ntot/chunksize
    if (ntot % chunksize) > 0:
        nchunk += 1

    print("nchunk:",nchunk)
    for i in xrange(nchunk):
        imin = i*chunksize
        imax = (i+1)*chunksize
        print("  %d %d (%d)" % (imin,imax,ntot))
        match = cols[zphot_matchname][imin:imax]
        
        nrows = match.size
        scinv = numpy.zeros((nrows, zlvals.size), dtype='f8') -9999.

        w,=numpy.where(match >= 0)
        print("    ",w.size,"with matches")
        if w.size > 0:
            # read the p(z) for the matches
            pofz = pzcols['pofz'][match[w]]

            pofz = pofz[:]*corr[:]

            for j in xrange(w.size):
                scinv[w[j], :] = scalc.calc_mean_scinv(pofz[j])
        # metadata only gets written once
        cols.write_column(colname, scinv, meta=meta)

def create_condor(type, procrun):
    import pbs
    import sdsspy
    proctype='regauss'
    procshort='rg'
    minscore=0.1

    condor_dir=path_join('~/condor/sweep_reduce',proctype,procrun)
    condor_dir = expand_path(condor_dir)
    if not os.path.exists(condor_dir):
        os.makedirs(condor_dir)

    win = sdsspy.window.Window()
    runs, reruns = win.runlist(minscore)

    runs.sort()
    print("writing to",condor_dir)
    for run in runs:
        for camcol in [1,2,3,4,5,6]:
            # create the condor file name
            rstr=sdsspy.files.run2string(run)
            job_name = '%s%s-%s-%s' % (procshort,type,run,camcol)

            script_base='%s%s-%s-%s-%s' % (procshort, type, procrun, rstr, camcol)
            condor_filename=script_base+'.condor'
            script_filename=script_base+'.py'
            condor_filename=path_join(condor_dir, condor_filename)
            script_filename=path_join(condor_dir, script_filename)

            condor_script="""
Universe        = vanilla
Notification    = Error
GetEnv          = True
Notify_user     = esheldon@bnl.gov
Requirements    = (CPU_Experiment == "astro") && (TotalSlots == 12 || TotalSlots == 8)
+Experiment     = "astro"
Initialdir      = {proc_dir}

Executable      = {script_base}.py
Output          = {script_base}.out
Error           = {script_base}.err
Log             = {script_base}.log

Queue\n""".format(proc_dir=condor_dir,
                  script_base=script_base)

            script="""#!{executable} -u
import os
os.environ['PHOTO_SWEEP'] = 'hdfs:///user/esheldon/boss/sweeps/dr8_final'
os.environ['SWEEP_REDUCE'] = 'hdfs:///user/esheldon/sweep-reduce'
from lensing.regauss import RegaussSweep
from es_sdsspy import sweeps
p=sweeps.Proc('{proctype}','{procrun}','{sweep_type}')
p.process_runs(RegaussSweep, runs={run}, camcols={camcol})
\n""".format(executable=sys.executable,
             proctype=proctype,
             procrun=procrun,
             sweep_type=type,
             run=run,
             camcol=camcol)


            print(condor_filename)
            with open(condor_filename,'w') as fobj:
                fobj.write(condor_script)
            with open(script_filename,'w') as fobj:
                fobj.write(script)
            os.system('chmod 755 '+script_filename)


def create_pbs(type, procrun):
    import pbs
    import sdsspy

    proctype='regauss'
    procshort='rg'
    minscore=0.1
    #queue = 'slow'
    #queue = 'batch'
    queue = 'fast'

    pbsdir=path_join('~/pbs',proctype,procrun)
    pbsdir = expand_path(pbsdir)
    if not os.path.exists(pbsdir):
        os.makedirs(pbsdir)

    espy_v = "-r ~esheldon/exports/espy-work"
    esutil_v="-r ~esheldon/exports/esutil-work"
    stomp_v="-r ~esheldon/exports/stomp-work"
    admom_v="-r ~esheldon/exports/admom-work"
    fimage_v="-r ~esheldon/exports/fimage-work"
    sdsspy_v="-r ~esheldon/exports/sdsspy-work"

    photo_redux   = os.environ['PHOTO_REDUX']
    photo_sweep   = os.environ["PHOTO_SWEEP"]
    photo_resolve = os.environ["PHOTO_RESOLVE"]
    photo_calib   = os.environ["PHOTO_CALIB"]

    setups = [espy_v,esutil_v,stomp_v,admom_v,fimage_v,sdsspy_v]
    setups = ['setup '+s for s in setups]

    setups += ['export PHOTO_REDUX=%s' % photo_redux,
               'export PHOTO_SWEEP=%s' % photo_sweep,
               'export PHOTO_CALIB=%s' % photo_calib,
               'export PHOTO_RESOLVE=%s' % photo_resolve]

    win = sdsspy.window.Window()
    runs, reruns = win.runlist(minscore)

    runs.sort()
    print("writing to",pbsdir)
    for run in runs:
        for camcol in [1,2,3,4,5,6]:
            # create the pbs file name
            rstr=sdsspy.files.run2string(run)
            job_name = '%s%s-%s-%s' % (procshort,type,run,camcol)

            pbsfilename='%s%s-%s-%s-%s.pbs' % (procshort, type, procrun, rstr, camcol)
            pbsfilename=path_join(pbsdir, pbsfilename)

            # the python commands to execute
            python_commands = """
from lensing.regauss import RegaussSweep
from es_sdsspy import sweeps
p=sweeps.Proc('%s','%s','%s')
p.process_runs(RegaussSweep, runs=%s, camcols=%s)""" % (proctype,procrun,type,run,camcol)

            ppy = pbs.PBSPython(pbsfilename,
                                python_commands,
                                job_name=job_name, 
                                setups=setups, 
                                queue=queue)
            ppy.write()


class RegaussSweep:
    def __init__(self, objs, **keys):
        self.objs = objs
        self.ra = RegaussAtlas()

        self.rmax = keys.get('rmax',22.)

        self.verbosity = keys.get('verbosity',1)

    def select(self):
        """
        """
        import es_sdsspy

        #raise RuntimeError("use mangle mask")
        print("Selecting objects")
        s = es_sdsspy.select.Selector(self.objs)
        print("  getting resolve logic")
        resolve_logic = s.resolve_logic()

        print("  getting mask logic")
        mask_logic = s.mask_logic('basic')

        print("  getting flag logic")
        flag_logic = s.flag_logic()

        print("  getting rmag logic")
        if 'devflux' in self.objs.dtype.names:
            rmag_logic = s.cmodelmag_logic("r", self.rmax)
        else:
            rmag_logic = s.modelmag_logic("r", self.rmax)

        logic = \
            resolve_logic & mask_logic & flag_logic & rmag_logic
        #logic = \
        #    resolve_logic & flag_logic & rmag_logic

        keep = where1(logic)
        print("  keeping %i/%i" % (keep.size, self.objs.size))

        if keep.size > 0:
            # original objs will live on with the caller
            self.objs = self.objs[keep]
        else:
            self.objs = None

    def process(self):
        import sdsspy

        self.select()
        if self.objs == None:
            return None

        ra = self.ra

        run=self.objs['run'][0]
        camcol=self.objs['camcol'][0]

        objs = self.objs
        self.make_output()

        b=eu.stat.Binner(objs['field'])
        b.dohist(binsize=1, rev=True)
        
        rev = b['rev']
        print("Processing objects by field")
        nf = b['hist'].size
        for i in xrange(nf):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                field = objs['field'][w[0]]
                print("  Processing %4i in %06i-%i-%04i  %3i/%3i" % (w.size,run,camcol,field,i+1,nf))
                
                for index in w:

                    obj = objs[index]
                    id = obj['id']
                    if ra.has_atlas(run,camcol,field,id):

                        if self.verbosity > 1:
                            print("    Processing %06i-%i-%04i-%04i" % (run,camcol,field,id))
                        self.output['has_atlas'][index] = 1

                        for fnum in [0,1,2,3,4]:
                            if self.verbosity > 1:
                                print("        ",sdsspy.FILTERCHAR[fnum])

                            rg = ra.regauss_obj(obj, fnum)
                            if rg is not None:
                                self.copy2output(index, rg, fnum)
        return self.output

    def copy2output(self, index, rg, fnum):
        
        # get a reference to this row, modifications to this will
        # also affect self.output
        data = self.output[index]

        if 'imstats' in rg:
            s = rg['imstats']
            if s is not None:
                data['wrow'][fnum]    = s['wrow']
                data['wcol'][fnum]    = s['wcol']
                data['Irr'][fnum]     = s['Irr']
                data['Irc'][fnum]     = s['Irc']
                data['Icc'][fnum]     = s['Icc']
                data['e1'][fnum]     = s['e1'] # just added these
                data['e2'][fnum]     = s['e2']
                data['a4'][fnum]      = s['a4']
                data['uncer'][fnum]   = s['uncer']
                data['s2n'][fnum]   = s['s2n']
                data['amflags'][fnum] = s['whyflag']
                data['amflags_str'][fnum] = s['whystr']
                data['numiter'][fnum] = s['numiter']
                
        if 'psfstats' in rg:
            s = rg['psfstats']
            if s is not None:
                data['Irr_psf'][fnum]     = s['Irr']
                data['Irc_psf'][fnum]     = s['Irc']
                data['Icc_psf'][fnum]     = s['Icc']
                data['a4_psf'][fnum]      = s['a4']
                data['amflags_psf'][fnum] = s['whyflag']
                data['amflags_psf_str'][fnum] = s['whystr']
 
        if 'corrstats' in rg:
            s = rg['corrstats']
            if s is not None:
                data['e1_lin'][fnum]        = s['e1']
                data['e2_lin'][fnum]        = s['e2']
                data['R_lin'][fnum]         = s['R']
                data['corrflags_lin'][fnum] = s['flags']


        if 'rgstats' in rg:
            s = rg['rgstats']
            if s is not None:
                data['f0flags'][fnum]    = s['f0flags']
                if s['f0flags'] == 0:
                    # we only get here if we were able to make f0
                    data['Irr_rg'][fnum]     = s['Irr']
                    data['Irc_rg'][fnum]     = s['Irc']
                    data['Icc_rg'][fnum]     = s['Icc']
                    data['a4_rg'][fnum]      = s['a4']
                    data['uncer_rg'][fnum]   = s['uncer']
                    data['amflags_rg'][fnum] = s['whyflag']
                    data['amflags_rg_str'][fnum] = s['whystr']
                    data['numiter_rg'][fnum] = s['numiter']
     

        if 'rgcorrstats' in rg:
            s = rg['rgcorrstats']
            if s is not None:
                data['e1_rg'][fnum]        = s['e1']
                data['e2_rg'][fnum]        = s['e2']
                data['R_rg'][fnum]         = s['R']
                data['corrflags_rg'][fnum] = s['flags']



    def make_output(self):
        import sdsspy

        print("creating output")
        dtype = self.output_dtype()
        self.output = numpy.zeros(self.objs.size, dtype=dtype)

        print("    copying from objs")
        eu.numpy_util.copy_fields(self.objs, self.output)

        mm=sdsspy.nmgy2mag(self.objs['modelflux'][:,2])
        self.output['modelmag_dered_r'] = mm-self.objs['extinction'][:,2]

        # set defaults
        output = self.output

        print("    setting defaults")
        # these get set to 9999
        flist = ['uncer','uncer_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = 9999

        # these get set to 2**15 = -32768 which means not processed
        flist = ['amflags','amflags_psf','amflags_rg',
                 'corrflags_lin','corrflags_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = 2**15

        for fnum in [0,1,2,3,4]:
            output['amflags_str'][:,fnum] = 'norun'
            output['amflags_psf_str'][:,fnum] = 'norun'
            output['amflags_rg_str'][:,fnum] = 'norun'

        # these get set to -9999
        flist = ['wrow','wcol',
                 'Irr','Irc','Icc','a4',
                 'Irr_psf','Irc_psf','Icc_psf','a4_psf',
                 'Irr_rg','Irc_rg','Icc_rg','a4_rg',
                 'e1_lin','e2_lin','R_lin',
                 'e1_rg','e2_rg','R_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = -9999

    def output_dtype(self):
        dt=[('run','i2'),
            ('rerun','i2'),
            ('camcol','i1'),
            ('field','i2'),
            ('id','i2'),
            ('thing_id','i4'),
            ('ra','f8'),  # do we need these?
            ('dec','f8'),
            ('modelmag_dered_r','f4'),

            ('has_atlas','i1'), # zero for no atlas

            ('wrow','5f8'),
            ('wcol','5f8'),
            ('Irr','5f8'),
            ('Irc','5f8'),
            ('Icc','5f8'),
            ('eq','5f8'),
            ('e1','5f8'), # just added these uncorrected shapes
            ('e2','5f8'),
            ('a4','5f8'),
            ('uncer','5f8'),
            ('s2n','5f8'),
            ('numiter','5i1'),
            ('amflags','5i2'),
            ('amflags_str','5S5'),

            ('Irr_psf','5f8'),
            ('Irc_psf','5f8'),
            ('Icc_psf','5f8'),
            ('a4_psf','5f8'),
            ('amflags_psf','5i2'),
            ('amflags_psf_str','5S5'),

            ('e1_lin','5f8'),
            ('e2_lin','5f8'),
            ('R_lin','5f8'),
            ('corrflags_lin','5i2'),  # flags during compea4 correction

            ('f0flags','5i1'),
            ('Irr_rg','5f8'),
            ('Irc_rg','5f8'),
            ('Icc_rg','5f8'),
            ('a4_rg','5f8'),
            ('uncer_rg','5f8'),
            ('numiter_rg','5i1'),
            ('amflags_rg','5i2'), # flags running admom on f0-epsilon in the rg process
            ('amflags_rg_str','5S5'),

            ('e1_rg','5f8'),    # corrected ellipticity
            ('e2_rg','5f8'),
            ('R_rg','5f8'),     # this is more stable than R_lin
            ('corrflags_rg','5i2')]  # flags during compea4 correction

        return dt


class RegaussAtlas:
    def __init__(self, verbose=False):
        self.verbose=verbose

        self.sweep_cache = SweepCache(verbose=self.verbose)
        self.atls_key = None
        self.field_key = None


    def regauss_id(self, type, run, camcol, field, id, filter):
        """
        Run regauss on a single object in a single band
        """

        obj = self.readobj(type, run, camcol, field, id)
        return self.regauss_obj(obj, filter)

    def regauss_obj(self, obj, filter):
        import admom
        from sdsspy.atlas.atlas import NoAtlasImageError

        run = obj['run']
        camcol = obj['camcol']
        field = obj['field']
        id = obj['id']

        c = sdsspy.FILTERNUM[filter]


        self.cache_psf(run, camcol, field)

        try:
            self.cache_atlas(run, camcol, field, id)
        except NoAtlasImageError:
            # we can safely ignore such errors, but we return
            # nothing
            return None

        im = numpy.array( self.atls['images'][c], dtype='f8')
        im -= self.atls['SOFT_BIAS']


        rowc=obj['rowc'][c]
        colc=obj['colc'][c]
        row=rowc - self.atls['row0'][c] - 0.5
        col=colc - self.atls['col0'][c] - 0.5

        sigsky = self.psfield['skysig'][0,c]
        psf = self.psfKLs[c].rec(rowc, colc, trim=True)

        # sometimes we get very negative pixels, not sure why
        #im=im.clip(-3*sigsky, im.max())

        guess_psf=admom.fwhm2mom( obj['psf_fwhm'][c], pixscale=0.4)/2.
        guess = guess_psf
        if 'm_rr_cc' in obj.dtype.names:
            if obj['m_rr_cc'][c] > guess_psf:
                guess=obj['m_rr_cc'][c]/2

        rg = admom.ReGauss(im, row, col, psf, 
                           sigsky=sigsky, guess=guess, guess_psf=guess_psf,
                           verbose=self.verbose)
        rg.do_all() 

        # keep these things in case we want to call show()
        self.rg=rg
        self.im=im
        self.psf=psf

        return rg

    def show(self, min=0, **keys):
        import fimage
        import biggles
        plt=images.view(self.im, show=False, min=min, **keys)
        s=self.rg['imstats']

        if s['whyflag'] == 0:
            levels=7

            wrow = s['wrow']
            wcol = s['wcol']

            model = fimage.model_image('gauss',
                                       self.im.shape,
                                       [wrow,wcol],
                                       [s['Irr'],s['Irc'],s['Icc']],
                                       counts=self.im.sum())

            cmod = biggles.Contours(model.transpose(), color='grey')
            cmod.levels = levels
            plt.add(cmod)
        plt.show()


    def readobj(self, type, run, camcol, field, id):
        obj = self.sweep_cache.readobj(type, run, camcol, field, id)
        return obj

    def has_atlas(self, run, camcol, field, id):
        from sdsspy.atlas.atlas import NoAtlasImageError

        if run == 1462 and camcol == 3 and field == 209:
            # this atlas image is corrupt
            return False
        try:
            self.cache_atlas(run, camcol, field, id)
            return True
        except NoAtlasImageError:
            # we can safely ignore such errors
            return False


    def cache_atlas(self, run, camcol, field, id):
        """
        This will cache reads so it can be used for multiple
        filter processings
        """
        import sdsspy
        key = '%06i-%i-%04i-%05i' % (run,camcol,field,id)
        if self.atls_key != key:
            self.atls = sdsspy.read('fpAtlas',run,camcol,field,id, 
                                    trim=True,
                                    verbose=self.verbose)
            self.atls_key = key

    def cache_psf(self, run, camcol, field):
        """

        This will cache reads so it can be used for multiple filter
        processings.  So process all objects from a field in order to be most
        efficient.

        """
        import sdsspy

        key = '%06i-%i-%04i' % (run,camcol,field)
        if self.field_key != key:
            self.psfield = sdsspy.read('psfield',run,camcol,field, 
                                       verbose=self.verbose)
            self.psfKLs = []
            for filter in sdsspy.FILTERCHARS:
                kl = sdsspy.read('psField',run,camcol,field,filter=filter,
                                 verbose=self.verbose)
                self.psfKLs.append(kl)

            self.field_key = key

def regauss_atlas(type, run, camcol, field, id, filter, 
                  show=False, showmin=0, verbose=False, **keys):
    """

    Run regauss on a given atlas image.  This is not optimal, since only a
    single flter is processed at a time, but is more for interactive work

    """

    rga = RegaussAtlas(verbose=verbose)
    rga.regauss1(type, run, camcol, field, id, filter)

    if show:
        rga.show(min=showmin)

    return rga.rg


def admom_atlas(type, run, camcol, field, id, filter, 
                objs=None,
                show=False, showmin=0, **keys):
    """
    Objs must be for this camcol!
    """
    import admom
    import sdsspy

    c = sdsspy.FILTERNUM[filter]

    if objs is None:
        sc=SweepCache()
        obj = sc.readobj(type, run, camcol, field, id)
    else:
        w=where1( (objs['field'] == field) & (objs['id'] == id) )
        if w.size == 0:
            raise ValueError("Object not found objs: "
                             "%06i-%i-%04i" % (run,camcol,id))
        obj = objs[w[0]]

    atls = sdsspy.read('fpAtlas',run,camcol,field,id, trim=True)
    im = numpy.array( atls['images'][c], dtype='f8') - atls['SOFT_BIAS']

    psfield = sdsspy.read('psField', run, camcol, field, lower=True)

    rowc=obj['rowc'][c]
    colc=obj['colc'][c]

    row=rowc - atls['row0'][c] - 0.5
    col=colc - atls['col0'][c] - 0.5

    sigsky = psfield['skysig'][0,c]
    Tguess_psf=admom.fwhm2mom( obj['psf_fwhm'][c], pixscale=0.4)
    if obj['m_rr_cc'][c] > Tguess_psf:
        Tguess=obj['m_rr_cc'][c]
    else:
        Tguess=Tguess_psf

    out = admom.admom(im,
                      row, 
                      col, 
                      sky=0.0,
                      sigsky=sigsky,
                      Tguess=Tguess)

    if show:
        import fimage
        import biggles
        plt=images.view(im, min=showmin, show=False, **keys)

        if out['whyflag'] == 0:
            levels=7

            wrow = out['wrow']
            wcol = out['wcol']

            model = fimage.model_image('gauss',
                                       im.shape,
                                       [wrow,wcol],
                                       [out['Irr'],out['Irc'],out['Icc']],
                                       counts=im.sum())


            cmod = biggles.Contours(model.transpose(), color='grey')
            cmod.levels = levels
            plt.add(cmod)
        plt.show()

    return out



class SweepCache:
    def __init__(self, type=None, data=None, verbose=False):
        self.verbose=verbose
        self.init(type=type, data=data)

    def init(self, type=None, data=None):
        if type is not None and data is not None:
            run=data['run'][0]
            camcol=data['camcol'][0]
            key = '%s-%06i-%d' % (type,run,camcol)
            self.key = key
            self.data = data
        else:
            self.key=None
            self.data=None


    def read(self,type,run,camcol, **keys):
        self.cache_column(type, run, camcol, **keys)
        return self.data

    def readfield(self, type, run, camcol, field, **keys):    
        self.cache_column(type, run, camcol, **keys)
        data = self.data
        w=where1(data['field'] == field)
        if w.size == 0:
            raise ValueError("field not found: %06i-%i-%04i" % (run,camcol,field))
        return data[w]

    def readobj(self, type, run, camcol, field, id, **keys):    
        self.cache_column(type, run, camcol, **keys)
        data = self.data
        w=where1( (data['field'] == field) & (data['id'] == id) )
        if w.size == 0:
            raise ValueError("object not found: %06i-%i-%04i-%05i" % (run,camcol,field,id))
        return data[w[0]]


    def cache_column(self, type,run,camcol, **keys):
        import sdsspy

        key = '%s-%06i-%d' % (type,run,camcol)
        if self.key != key:
            data = sdsspy.read('calibobj.%s' % type,
                               run,camcol, lower=True, verbose=self.verbose, 
                               **keys)
            self.key = key
            self.data = data


class Collator:
    """

    Collate the regauss outputs into columns.

    Note we put both types ('gal','star') into the same run.  This means if you
    make indices you'll have to re-make them when the other type is added.

    """
    
    def __init__(self, procrun, sweeptype, fs='nfs', coldir=None):
        import es_sdsspy

        self.proctype='regauss'
        self.sweeptype=sweeptype
        self.procrun=procrun
        self.fs = fs
        self.coldir=coldir

        self.sweep_cols = es_sdsspy.sweeps.open_columns(sweeptype)



    def collate_as_columns_byband(self):
        """
        Collate the collated files from the run into a columns directory.
        Multi-band columns are split to, e.g., flux_r, flux_i, ...
        """

        c = self.open_columns()

        print("Will write in coldir:", c.dir)
        if c.dir_exists():
            raise ValueError("Columns already exist, start from scratch")
        c.create()

        files = self.file_list()

        # create the meta column from input data
        tmp,header=eu.io.read(files[0],header=True)
        meta = self.make_meta(header)
        c.write_column('meta', meta)

        for f in files:
            print("Reading:",f)
            
            tmp=eu.io.read(f)

            if len(tmp) == 0:
                print("No results")
            else:
                print("  Found",tmp.size,"... Creating output")
                out_dict = self.create_output(tmp)

                print("  Writing columns")
                for name in out_dict:
                    c.write_column(name, out_dict[name])
        print("Done")



    def create_output(self, st):
        import sdsspy

        bands = ['u','g','r','i','z']

        out_dict = {}
        for d in st.dtype.descr:
            name = str( d[0] )

            if len(d) == 3:
                if isinstance(d[2],tuple):
                    sz=d[2][0]
                else:
                    sz=d[2]
                if sz == 5:
                    for bandi in xrange(5):
                        fname = name+'_'+bands[bandi]
                        out_dict[fname] = st[name][:, bandi]
            else:
                out_dict[name] = st[name]

        # match up to sweeps and get some info
        print("    matching to primary sweeps")
        w=self.sweep_cols['thing_id'].match(st['thing_id'])
        #print("    found %s/%s" % (w.size,st.size))
        if w.size != st.size:
            raise ValueError("Did not match all")

        out_dict['primgal_id'] = w

        print("    Getting psf ellip")
        Irr_psf = st['Irr_psf']
        Irc_psf = st['Irc_psf']
        Icc_psf = st['Icc_psf']
        T_psf = st['Irr_psf'] + st['Icc_psf']

        if 'e1' not in st.dtype.names:
            e1,e2 = self.make_e1e2(st)

        print("    Copying ellip,mags")
        for f in sdsspy.FILTERCHARS:
            fnum = sdsspy.FILTERNUM[f]
            out_dict['e1_psf_'+f] = (Icc_psf[:,fnum]-Irr_psf[:,fnum])/T_psf[:,fnum]
            out_dict['e2_psf_'+f] = 2*Irc_psf[:,fnum]/T_psf[:,fnum]

            out_dict['modelmag_dered_'+f] = self.sweep_cols['modelmag_dered_'+f][w]
            ext = self.sweep_cols['extinction_'+f][w]
            out_dict['extinction_'+f] = ext

            cmodelname = 'cmodelmag_dered_'+f
            if cmodelname in self.sweep_cols:
                out_dict[cmodelname] = self.sweep_cols[cmodelname][w]

            if self.sweeptype == 'star':
                flux = self.sweep_cols['psfflux_'+f][w]
                ivar = self.sweep_cols['psfflux_ivar_'+f][w]
                flux_dered = sdsspy.dered_fluxes(ext, flux)
                mag_dered = sdsspy.nmgy2mag(flux_dered)

                out_dict['psfflux_'+f] = flux
                out_dict['psfflux_ivar_'+f] = ivar
                out_dict['psfmag_dered_'+f] = mag_dered

            if 'e1' not in st.dtype.names:
                out_dict['e1_'+f] = e1[:,fnum]
                out_dict['e2_'+f] = e2[:,fnum]


        out_dict['photoid'] = sdsspy.photoid(st)

        return out_dict

    def replace_e1e2(self):
        """
        I had a bug calculating e1,e2 first time.  This will replace it
        """
        c = self.open_columns()
        for band in ['u','g','r','i','z']:
            print('band:',band)
            e1col = 'e1_'+band
            e2col = 'e2_'+band

            print("  reading Irr")
            Irr = c['Irr_'+band][:]
            print("  reading Irc")
            Irc = c['Irc_'+band][:]
            print("  reading Icc")
            Icc = c['Icc_'+band][:]
            print("  reading flags")
            flags = c['amflags_'+band][:]

            T = Irr + Icc

            e1 = 0*Irr.copy() -9999
            e2 = e1.copy()

            w,=numpy.where( (flags==0) & (T > 0) )
            if w.size > 0:
                e1[w] = (Icc[w]-Irr[w])/T[w]
                e2[w] = 2*Irc[w]/T[w]

            print("  writing new column:",e1col)
            c.write_column(e1col, e1, create=True)
            print("  writing new column:",e2col)
            c.write_column(e2col, e2, create=True)

    def make_e1e2(self, st):
        """
        We forgot to put regular e1,e2 in struct
        """
        Irr = st['Irr']
        Irc = st['Irc']
        Icc = st['Icc']
        T = Irr + Icc

        e1 = 0*Irr.copy() -9999
        e2 = e1.copy()

        wg=numpy.where( T > 0 )
        if wg[0].size > 0:
            e1[wg] = (Icc[wg]-Irr[wg])/T[wg]
            e2[wg] = 2*Irc[wg]/T[wg]
        return e1,e2


    def add_rotated_e1e2(self, filters=['u','g','r','i','z'], 
                         system='eq', 
                         detrend=False,
                         rmag_max=21.8):
        from . import rotation
        rotator = rotation.SDSSRotator(system)
        c=self.open_columns()
        print("reading id info")
        runs    = c['run'][:]
        camcols = c['camcol'][:]
        fields  = c['field'][:]

        if self.sweeptype == 'gal':
            if detrend:
                rmstr='%0.1f' % rmag_max
                rmstr = rmstr.replace('.','')
                e1name = 'e1_rg_dt'+rmstr
                e2name = 'e2_rg_dt'+rmstr
                flagname = 'dt'+rmstr+'_flag'
                rotname='rot_dt'+rmstr
            else:
                e1name = 'e1_rg'
                e2name = 'e2_rg'
                flagname = 'corrflags_rg'
                rotname='rot'
        else:
            e1name = 'e1'
            e2name = 'e2'
            flagname = 'amflags'

        for filter in filters:
            print("filter: '%s'" % filter)

            e1col=e1name + '_'+system+'_'+filter
            e2col=e2name + '_'+system+'_'+filter
            rotcol=rotname+'_'+system+'_'+filter

            print("  Reading flags,e1,e2")

            flags = c[flagname+'_'+filter][:]
            e1pix = c[e1name+'_'+filter][:]
            e2pix = c[e2name+'_'+filter][:]
            print("  selecting flags == 0")
            w=where1(flags == 0)
            print("    found: %i/%i" % (w.size,flags.size))
            print("  rotating")
            te1, te2, tangle = rotator.rotate(runs[w], camcols[w], fields[w], filter,
                                              e1pix[w], e2pix[w], getrot=True)

            e1new = numpy.zeros(e1pix.size, dtype='f8') - 9999.
            e2new = e1new.copy()
            angles = e1new.copy()

            e1new[w] = te1
            e2new[w] = te2
            angles[w] = tangle

            print("writing new column:",e1col)
            c.write_column(e1col, e1new, create=True)
            print("writing new column:",e2col)
            c.write_column(e2col, e2new, create=True)
            print("writing rotation column:",rotcol)
            c.write_column(rotcol, angles, create=True)


    def create_indices(self):
        c = self.open_columns()
        print('creating indices in coldir:',c.dir)

        if self.sweeptype == 'gal':
            cnames = ['run','camcol','thing_id','photoid',
                      'cmodelmag_dered_r',
                      'amflags_r','amflags_i',
                      'amflags_rg_r','amflags_rg_i',
                      'corrflags_lin',
                      'corrflags_rg',
                      'R_lin_r','R_lin_i',
                      'R_rg_r','R_rg_i',
                      'e1_lin_r','e1_lin_i',
                      'e2_lin_r','e2_lin_i',
                      'e1_rg_r','e1_rg_i',
                      'e2_rg_r','e2_rg_i']
        else:
            cnames = ['run','camcol','thing_id','photoid',
                      'modelmag_dered_r',
                      'psfmag_dered_r',
                      'amflags_r','amflags_i',
                      'e1_r','e2_r',
                      'e1_i','e2_i']

        for n in cnames:
            if n in c:
                print("Creating index for:",n)
                c[n].create_index(force=True)



    def make_meta(self, header):
        import datetime
        date = str( datetime.datetime.now() )
        out = {}
        tags=['photoop_v',
              'idlutils_v',
              'photo_sweep',
              'photo_resolve',
              'photo_calib']
        print("ADD OTHER SOFTWARE VERSIONS TO HEADER")
        for n in tags:
            if n in header:
                out[n] = str(header[n])
        out['date'] = date
        out['proctype'] = self.proctype
        out['procrun'] = self.procrun
        out['sweeptype'] = self.sweeptype

  
        return out


    



    def open_columns(self):
        import columns
        d = self.columns_dir()
        return columns.Columns(d)

    def columns_dir(self):
        import es_sdsspy

        if self.coldir is None:
            p=es_sdsspy.sweeps.Proc(self.proctype,self.procrun,self.sweeptype) 
            d = p.columns_dir()
        else:
            d=self.coldir
        return d

    def output_dir(self):
        import es_sdsspy

        if self.fs == 'hdfs':
            d = 'hdfs:///user/esheldon/sweep-reduce/regauss/%s' % self.procrun
        else:
            p=es_sdsspy.sweeps.Proc(self.proctype,self.procrun,self.sweeptype) 
            d = p.output_dir()
        return d

    def file_list(self):
        indir = self.output_dir()
        pattern='regauss-{sweeptype}-{procrun}-*-*-*.rec'
        pattern=pattern.format(sweeptype=self.sweeptype, procrun=self.procrun)

        pattern = os.path.join(indir,pattern)
        print("Searching for pattern:",pattern)
        if pattern[0:4] == 'hdfs':
            files=eu.hdfs.ls(pattern)
            files = ['hdfs://' + f for f in files]
        else:
            files = glob.glob(pattern)

        files.sort()
        return files





def add_radec_to_columns(procrun):
    import es_sdsspy

    gals = es_sdsspy.sweeps.open_columns('gal')
    c = open_columns(procrun)

    if 'ra' in c or 'dec' in c:
        raise ValueError("ra or dec already in columns: "+c.dir)

    print("Reading primgal_id for correspondence")
    galid = c['primgal_id'][:]

    print("Reading ra")
    ra = gals['ra'][galid]
    print("Writing matching ra")
    c.write_column('ra', ra)

    print("Reading dec")
    dec = gals['dec'][galid]
    print("Writing matching dec")
    c.write_column('dec', dec)

def match_yale(procrun):
    import es_sdsspy

    dir = '~/lensing/yale'
    yfile=path_join(dir,'sdss_images1_iv.rec')
    outfile = path_join(dir, 'match-yale-%s.rec' % procrun)
    print("Will write to output file:",outfile)

    yale = eu.io.read(yfile, verbose=True)

    urun = numpy.unique(yale['run'])

    max_rmag=22

    c=open_columns(procrun)
    append=False
    for run in urun:
        
        wyale = where1(yale['run'] == run)

        print('-'*70)
        print("Matching run:",run)
        wrun = c['run'] == run

        print("Found",wrun.size,"in run")

        print("Getting whyflag== 0 and r <",max_rmag)
        whyflag = c['whyflag_rg_r'][wrun]
        mag = c['cmodelmag_dered_r'][wrun]

        w=where1( (whyflag == 0) & (mag < max_rmag) )

        w = wrun[w]

        print("Found:",w.size)

        print("Reading ra")
        ra=c['ra'][w]
        print("Reading dec")
        dec=c['dec'][w]

        h = eu.htm.HTM()

        matchrad = 1.0
        print("Matching ra/dec within",matchrad,"arcsec")
        myale, mrg, dist = h.match(yale['ra'][wyale], yale['dec'][wyale], ra, dec, matchrad/1000.)

        if myale.size == 0:
            raise ValueError("no matches!")
        print("Found %d/%d ra/dec matches" % (mrg.size, wyale.size))

        # index back into the original columns
        w = w[mrg]

        columns = ['run','rerun','camcol','field','id',
                   'ra','dec','e1_rg_r','e2_rg_r','momerr_r','r_r',
                   'cmodelmag_dered_r']
        data = c.read_columns(columns, rows=w)

        gc = es_sdsspy.sweeps.open_columns('gal')
        meta = gc['meta'].read()

        hdr={'regauss_procrun':procrun,
             'photo_sweep': meta['photo_sweep'],
             'matchrad':matchrad,
             'max_rmag':max_rmag}

        eu.io.write(outfile, data, header=hdr, append=append, delim=' ')
        append=True
