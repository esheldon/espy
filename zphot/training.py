"""
See docs in weighting.py
"""
from __future__ import print_function
import os
from sys import stdout,stderr
import pprint
import copy

import numpy
from numpy import where, sqrt, linspace, exp, pi as PI

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

import sdsspy

import biggles

import zphot

class Primus:
    """
    Special class to trim down the primus sample
    """
    def __init__(self):
        pzdir=zphot.photoz_dir()
        self.dir = path_join(pzdir,'training','original-2010-09-19')
        self.infile = path_join(self.dir,'primus.zerod.10oct29.fits.gz')
        self.outfile = path_join(self.dir,'primus.zerod.10oct29.zconf4.rec')

        self.data = eu.io.read(self.infile,lower=True,verbose=True)
        eu.numpy_util.to_native(self.data, inplace=True)

    def get_trimmed(self):
        w=where1(self.data['zprimus_zconf'] == 4)
        trimmed = self.data[w]
        return trimmed
    def write_trimmed(self):
        trimmed = self.get_trimmed()
        eu.io.write(self.outfile, trimmed, verbose=True)

    def plot_trimmed(self, show=False):
        import converter

        binsize = 0.02
        bmin = 0.01
        bmax = 1.25
        psfile = self.outfile.replace('.rec','-comparez.eps')

        plt = biggles.FramedPlot()
        plt.xlabel = 'z'

        hall = eu.stat.histogram(self.data['zprimus'], binsize=0.02, min=bmin, max=bmax, more=True)
        phall=biggles.Histogram(hall['hist'], x0=hall['low'][0], binsize=binsize)
        phall.label = 'all'

        trimmed = self.get_trimmed()
        htrim = eu.stat.histogram(trimmed['zprimus'], binsize=0.02, min=bmin, max=bmax, more=True)
        phtrim=biggles.Histogram(htrim['hist'], x0=htrim['low'][0], binsize=binsize, color='red')
        phtrim.label = 'zconf = 4'

        key = biggles.PlotKey(0.6,0.9,[phall,phtrim])

        plt.add(phall, phtrim, key)

        if show:
            plt.show()
        plt.write_eps(psfile)
        converter.convert(psfile,dpi=90,verbose=True)



class Training:
    """
    Code to match the training set against an imaging data set

Okay, I've checked all the files.
tkrs,2slaq,sdss, and cfrs are ready to go.
The others need the following cuts.
deep2: zqual >= 3
zcosmos: cc= 3.4 || 3.5 || 4.4. || 4.5 || 9.5
vvds: zqual=3 || 4

If you want to clean the files up before or after you send them to me,
doesn't matter. Before getting the p(z)'s I suggest esimating the weighted
N(z) for one of the magnitude limited samples, such as zcosmos and vvds, to
see if things are decent. I can do that part if you're not comfortable with
it.
Carlos

    """
    def __init__(self, train_sample=None, no_photo_cuts=False):
        self.train_sample = train_sample
        self.dir_matched = None
        self.no_photo_cuts = no_photo_cuts

        pzdir = zphot.photoz_dir()
        self.basedir = path_join(pzdir, 'training')
        self.prefix = 'train'

        self.conf=None

        if self.train_sample is not None:
            self.init(self.train_sample)

    def init(self, train_sample):
        self.types = {'primus':'primus.zerod.10oct29.zconf4',
                      '2slaq':'2slaq',
                      'cfrs':'cfrs',
                      'cnoc2':'cnoc2.cut',
                      'deep2':'deep2.form.fix',
                      'sdss':'sdssobjids',
                      'tkrs':'tkrs-fix',
                      'vvds':'vvds',
                      'zcosmos':'zcosmos'}

        self.dir_matched = path_join(self.basedir, 
                                     'matched', 
                                     train_sample)
        if self.no_photo_cuts:
            self.dir_matched = path_join(self.dir_matched,'no_photo_cuts')

        self.conf = zphot.read_config('train',self.train_sample)

        self.photo_conf = zphot.read_config('zinput',self.conf['photo_sample'])
        pprint.pprint(self.conf)

    def matched_dtype(self, dtype_in, type):
        dtype = numpy.dtype(dtype_in)

        dt = copy.deepcopy(dtype.descr)
        dt.append( ('photoid','i8') )

        dt.append( ('mspec','i4') )
        dt.append( ('mphot','i4') )

        dt.append( ('modelmag_dered_u', 'f4') )
        dt.append( ('modelmag_dered_g', 'f4') )
        dt.append( ('modelmag_dered_r', 'f4') )
        dt.append( ('modelmag_dered_i', 'f4') )
        dt.append( ('modelmag_dered_z', 'f4') )

        dt.append( ('cmodelmag_dered_r','f4') )

        '''
        dt.append( ('modelmag_dered_err_u', 'f4') )
        dt.append( ('modelmag_dered_err_g', 'f4') )
        dt.append( ('modelmag_dered_err_r', 'f4') )
        dt.append( ('modelmag_dered_err_i', 'f4') )
        dt.append( ('modelmag_dered_err_z', 'f4') )
        '''

        dt.append( ('survey_primary','i1') )
        dt.append( ('cmodelmag_dered_err_r','f4') )

        dt.append( ('psf_fwhm_r','f4') )

        if type == 'primus':
            dt.append( ('z','f4') )

        names = [d[0] for d in dt]

        self.extra_colnames = []
        if 'extra_columns' in self.photo_conf:
            for cdict in self.photo_conf['extra_columns']:
                if cdict['name'] not in names:
                    dt.append( (cdict['name'], cdict['dtype']) )
                    self.extra_colnames.append(cdict['name'])
                else:
                    print("extra columns '%s' already in data" % cdict['name'])

        return dt

    def match(self):
        if self.conf is None:
            raise ValueError("Init with sample name")

        if not os.path.exists(self.dir_matched):
            os.makedirs(self.dir_matched)

        # run the same selection as we did in making
        # the photoz inputs.  Note I'm *not* putting
        # the mask or rmax cuts here, that should come
        # out in the weighting

        photo_sample = self.conf["photo_sample"]
        zcs = zphot.select.ColumnSelector(photo_sample)

        zcs.select()

        indices = zcs.keep_indices.copy()

        print("\nReading photoid,mag,cmag")
        photoid = zcs.cols['photoid'][indices]

        mag_u    = zcs.cols['modelmag_dered_u'][indices]
        mag_g    = zcs.cols['modelmag_dered_g'][indices]
        mag_r    = zcs.cols['modelmag_dered_r'][indices]
        mag_i    = zcs.cols['modelmag_dered_i'][indices]
        mag_z    = zcs.cols['modelmag_dered_z'][indices]


        cmag_r = zcs.cols['cmodelmag_dered_r'][indices]
        cmag_err_r = zcs.cols['cmodelmag_dered_err_r'][indices]

        if 'survey_primary' in zcs.cols:
            # we are doing multi-epoch matching
            print("Reading survey_primary, psf_fwhm_r")
            survey_primary = zcs.cols['survey_primary'][indices]
        else:
            survey_primary = None

        psf_fwhm_r = zcs.cols['psf_fwhm_r'][indices]

        extra_columns = None

        rmin = zcs.conf['cmodel_rmin']
        rmax = zcs.conf['cmodel_rmax']

        if (cmag_r.min() < rmin) or (cmag_r.max() > rmax):
            raise ValueError("Expected max rmag within [%s,%s]\n" % (rmin,rmax))

        # get all ra/dec for matching
        print('Reading ra/dec')
        ra = zcs.cols['ra'][indices]
        dec = zcs.cols['dec'][indices]

        h = eu.htm.HTM(10)
        htmid = h.lookup_id(ra,dec)
        minid = htmid.min()
        maxid = htmid.max()
        print("Getting htmrev")
        hist,htmrev = eu.stat.histogram(htmid-minid, rev=True)

        htminfo={'htmid':htmid,
                 'minid':minid,
                 'maxid':maxid,
                 'htmrev':htmrev}

        # match within 2 arcsec
        #matchrad_arcsec = 2.0
        #matchrad_arcsec = 1.0
        #matchrad = matchrad_arcsec/3600.0
        
        for type in self.types:
            print('\n' + '-'*70)
            spec = self.read_original(type)

            w = self.make_specific_cuts(spec, type)
            spec = spec[w]

            mspec, mphot = self.extract_good_matches(ra, dec, cmag_r, cmag_err_r, survey_primary, htminfo,
                                                     spec['ra'], spec['dec'])

            mspec_u = numpy.unique(mspec)
            print('%s/%s of spec matches are unique' % (mspec_u.size,mspec.size))
            mphot_u = numpy.unique(mphot)
            print('%s/%s of photo matches are unique' % (mphot_u.size,mphot.size))

            dt = self.matched_dtype(spec.dtype.descr, type)

            # need call to matched_dtype to establish extra columns
            if extra_columns is None:
                extra_columns={}
                if len(self.extra_colnames) > 0:
                    for n in self.extra_colnames:
                        print("Reading extra:",n)
                        extra_columns[n] = zcs.cols[n][indices]


            output = numpy.zeros(mspec.size,dtype=dt)
            eu.numpy_util.copy_fields(spec[mspec],output)

            output['mspec'] = mspec
            output['mphot'] = mphot

            if 'primus' in type:
                output['z'] = output['zprimus']
            
            print("copying photoid")
            output['photoid'] = photoid[mphot]

            print("copying mags/fwhm")

            # note we don't want to read the matches directly from the columns
            # database this way because there are duplicates in mphot.  We must
            # keep in memory to do duplicate subscripting

            output['modelmag_dered_u']  = mag_u[mphot]
            output['modelmag_dered_g']  = mag_g[mphot]
            output['modelmag_dered_r']  = mag_r[mphot]
            output['modelmag_dered_i']  = mag_i[mphot]
            output['modelmag_dered_z']  = mag_z[mphot]

            output['cmodelmag_dered_r'] = cmag_r[mphot]
            output['cmodelmag_dered_err_r'] = cmag_err_r[mphot]

            output['psf_fwhm_r'] = psf_fwhm_r[mphot]

            if survey_primary is None:
                output['survey_primary'] = numpy.ones(mphot.size,dtype='i1')
            else:
                output['survey_primary'] = survey_primary[mphot]

            for n in extra_columns:
                output[n] = extra_columns[n][mphot]


            # sanity check
            if not self.no_photo_cuts:
                if ((output['cmodelmag_dered_r'].max() > rmax) or
                    (output['cmodelmag_dered_r'].min() < rmin) ):
                    raise ValueError("Expected output rmag in [%s,%s]\n" % (rmin,rmax))

            # during matching we might end up with some too-faint mags for
            # the sdss
            if type == 'sdssobjids':
                print("    limiting SDSS training set to r < %s" \
                             % self.conf['sdss_rmax'])
                w,=where( output['cmodelmag_dered_r'] < self.conf['sdss_rmax'])
                print("        keeping %s/%s" % (w.size,mphot.size))
                output = output[w]
            

            # this may be sample dependant.
            zmin = self.conf['zmin']
            zmax = self.conf['zmax']
            print("Selecting spec to be in redshift "
                  "range: [%s,%s]" % (zmin,zmax))
            z = output['z']
            w,=where((z > zmin) & (z < zmax))
            print("Keeping %s/%s" % (w.size,mphot.size))
            output = output[w]



            outname=self.fname_matched(type)
            print("Writing to match file:",outname)
            hdr = copy.deepcopy(self.conf)

            origfile=self.fname_original(type)
            hdr['original_file'] = os.path.basename(origfile)
            eu.io.write(outname, output, delim=' ', header=hdr)

    def extract_good_matches(self, ra, dec, cmag_r, cmag_err_r, survey_primary, htminfo,
                             spec_ra, spec_dec, 
                             nsigma=2.5, matchrad=1./3600.):
        """
        Match by ra/dec.

        Demand one of the matches is survey_primary, and use it as a reference
        point to remove outliers at nsigma in flux.

        Choosing nsigma=2.5 and 1'' match radius to make sure we get good matches

        """

        if survey_primary is None:
            maxmatch=1
        else:
            maxmatch=100
        print("Matching to %0.2f arcsec" % (matchrad*3600.,))
        h=eu.htm.HTM(10)

        mspec,mphot,d12 = h.match(spec_ra,spec_dec, ra, dec,
                                  matchrad,
                                  htmid2=htminfo['htmid'], 
                                  minid=htminfo['minid'],
                                  maxid=htminfo['maxid'],
                                  htmrev2=htminfo['htmrev'],
                                  maxmatch=maxmatch)

        print("Matched: %s/%s" % (mspec.size,spec_ra.size))
        if mspec.size == 0 or survey_primary==None:
            return mspec, mphot


        print("  unique:",numpy.unique(mspec).size)
        print("  getting good flux matches")

        h,rev = eu.stat.histogram(mspec, binsize=1, rev=True)
        keep = numpy.zeros(mspec.size, dtype='i1')

        for i in xrange(h.size):
            if rev[i] != rev[i+1]:

                # indices of mspec/mphot that have a particular index into
                # mspec
                w=rev[ rev[i]:rev[i+1] ]


                # make sure this makes sense
                wbad=where1(mspec[w] != mspec[w[0]])
                if wbad != 0:
                    raise ValueError("expected all indices to have same mspec")
                
                # indices into the photometric sample
                wphot=mphot[w]
                flux, ivar = sdsspy.mag2nmgy(cmag_r[wphot], cmag_err_r[wphot])

                werr=where1(ivar > 0.)
                if werr.size > 0:
                    # reduce the indices into mspec/mphot
                    w=w[werr]

                    # trim photometric data
                    flux=flux[werr]
                    ivar=ivar[werr]
                    wphot=wphot[werr]


                    # demand a survey primary match
                    wp = where1(survey_primary[wphot] == 1)
                    if wp.size > 0:
                        # always keep the primary
                        keep[w[wp]] = 1

                        wphot_primary = wphot[wp[0]]


                        # flux in primary will be reference
                        pflux, pivar = sdsspy.mag2nmgy(cmag_r[wphot_primary], 
                                                       cmag_err_r[wphot_primary])

                        # nsigma in the combined error
                        err = sqrt(1./ivar)
                        perr = sqrt(1./pivar)

                        comb_err = sqrt(err**2 + perr**2)

                        wgood = where1( abs(flux-pflux) < nsigma*comb_err)
                        if wgood.size > 0:
                            w=w[wgood]
                            keep[w] = 1

        wkeep = where1(keep == 1)
        print("  kept: %s/%s" % (wkeep.size, keep.size))
        if wkeep.size > 0:
            mspec = mspec[wkeep]
            mphot = mphot[wkeep]
        else:
            mspec = numpy.array([],dtype='i4')
            mphot = numpy.array([],dtype='i4')

        return mspec, mphot

    def plot_matches(self, type, index=None, data=None, nmatchmin=5):
        if data is None:
            data=self.read_matched(type)

        if index is None:
            index=numpy.arange(data.size, dtype='i4')
        else:
            index=numpy.array(index, ndmin=1, dtype='i4')

        mspec = data['mspec'][index]
        h,rev = eu.stat.histogram(mspec, binsize=1, rev=True)

        for i in xrange(h.size):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                if w.size >= nmatchmin:
                    cmag = data['cmodelmag_dered_r'][w]
                    cmagerr = data['cmodelmag_dered_err_r'][w]
                    flux, ivar = sdsspy.mag2nmgy(cmag, cmagerr)
                    err=sqrt(1./ivar)

                    minf = (flux-3.5*err).min()
                    maxf = (flux+3.5*err).max()
                    fvals = linspace(minf, maxf, 100)

                    plt=biggles.FramedPlot()

                    for wi in xrange(w.size):
                        if data['survey_primary'][w[wi]] == 1:
                            color='red'
                            pmag = cmag[wi]
                        else:
                            color='black'
                        g = exp(-0.5*( (fvals-flux[wi])/err[wi])**2)/sqrt(2.*PI)/err[wi]
                        gcurve = biggles.Curve(fvals, g, color=color)
                        plt.add(gcurve)

                    #std = cmag.std()
                    #binsize=std*0.4
                    #plt=eu.plotting.bhist(cmag, binsize=binsize, show=False)

                    nlab=biggles.PlotLabel(0.05,0.95,'nmatch: %s' % w.size, halign='left')
                    maglab = biggles.PlotLabel(0.05,0.9,'primary mag: %0.2f' % pmag, halign='left')
                    plt.add(nlab, maglab)
                    plt.title=type
                    plt.xlabel = 'flux'
                    plt.show()

                    k=raw_input('hit a key (q to quit): ')
                    if k.lower() == 'q':
                        return
                else:
                    print("less than %s matches" % nmatchmin)
            


    def make_specific_cuts(self, data, type):
        if type == 'deep2':
            w,=where(data['zqual'] >= 3)
            print("deep2: keeping %s/%s" % (w.size, data.size))
            return w
        elif type == 'zcosmos':
            cc = data['cc']
            w,=where((cc==3.4) | (cc==3.5) | (cc==4.4) | (cc==4.5) | (cc==9.5))
            print("zcosmos: keeping %s/%s" % (w.size, data.size))
            return w
        elif type == 'vvds':
            zflag = data['zflag']
            w,=where((zflag==3) | (zflag==4))
            print("vvds: keeping %s/%s" % (w.size, data.size))
            return w
        else:
            return numpy.arange(data.size,dtype='i4')



    def read_original(self, type):
        fname = self.fname_original(type)
        print('Reading original spec file:', fname)
        return eu.io.read(fname)

    def read_matched(self, type):
        fname = self.fname_matched(type)
        print('Reading matched spec file:', fname)
        return eu.io.read(fname)


    def fname_original(self, type):
        ext='.rec'
        dir = path_join(self.basedir, 'original-2010-09-19')
        long_type = self.types[type]
        return os.path.join(dir, long_type+ext)

    def fname_matched(self, type):
        dir = self.dir_matched
        if type == 'all':
            long_type='all'
        else:
            long_type = self.types[type]
        fname = long_type+'-match-'+self.train_sample+'.rec'
        fname = os.path.join(dir, fname)
        return fname

    def plotdir(self):
        dir = self.dir_matched
        return path_join(dir, 'plots')

    def seeing_plotfile(self, name):
        dir = self.plotdir()
        fname = '%s-match-seeing-%s.eps' % (name,self.train_sample)
        fname=path_join(dir, fname)
        return fname

    def plotfile(self, type, extra=None):
        dir = self.plotdir()
        fname = self.fname_matched(type)
        fname = os.path.basename(fname)
        fname = fname.replace('.rec','.eps')
        if extra is not None:
            fname = fname.replace('.eps','-'+extra+'.eps')
        fname = path_join(dir,fname)
        return fname

    def plot_all_radec(self):
        for t in self.types:
            self.plot_radec(t)

    def plot_radec(self, type):
        """
        ra/dec plot of all points and the matched points
        """
        import biggles
        import converter
        from biggles import FramedPlot,Points

        print()
        dir=self.plotdir()
        if not os.path.exists(dir):
            os.makedirs(dir)

        psfile = self.plotfile(type)

        orig = self.read_original(type)
        mat = self.read_matched(type)

        plt=FramedPlot() 

        symsize = 2
        if type == 'sdss' or type == 'other':
            symsize = 0.25

        plt.add( Points(orig['ra'],orig['dec'],type='dot', size=symsize) )

        plt.add( Points(mat['ra'],mat['dec'],type='dot',
                        color='red', size=symsize) )
        plt.xlabel = 'RA'
        plt.ylabel = 'DEC'
        print("Writing eps file:", psfile)
        plt.write_eps(psfile)
        converter.convert(psfile, dpi=120, verbose=True)


    def plot_seeing(self, types=['primus','vvds','zcosmos','deep2'], 
                    yrange=None):
        """

        The BOSS all is normalized to one
        The summed is normalized to one

        Others are normalized relative to the
        summed

        """
        import pcolors
        import es_sdsspy
        import converter

        binsize=0.025

        if types is None or types is 'all':
            name='all'
            # convert keys to a list
            types=list(self.types)
        else:
            name='-'.join(types)


        plt=biggles.FramedPlot()
        ntype=len(types)
        #colors=pcolors.rainbow(ntype, 'hex')
        colors=pcolors.rainbow(ntype+1, 'hex')

        allhist=[]
        
        c = es_sdsspy.sweeps_collate.open_columns('primgal')
        print("Reading all psf_fwhm")
        print("  checking window")
        inbasic = c['inbasic'][:]
        psf_fwhm = c['psf_fwhm_r'][:]
        mag = c['cmodelmag_dered_r'][:]
        print("  psf_fwhm_r shape:",psf_fwhm.shape)

        w=where1((inbasic == 1) & (mag < 21.8))
        print("  keeping: %s/%s" %(w.size,c['inbasic'].size))
        psf_fwhm = psf_fwhm[w]
        print("median seeing:",numpy.median(psf_fwhm))
        
        b=eu.stat.Binner(psf_fwhm)
        xmin=0.6
        xmax=1.8
        b.dohist(binsize=binsize, min=xmin, max=xmax)
        b.calc_stats()
        h = b['hist']/float(b['hist'].sum())
        #h = 0.3*b['hist']/float(b['hist'].max())
        x = b['center']
        htot = h.copy()
        htot[:] = 0

        print("  histogramming")
        # regular histogram just for the key
        ph = biggles.Histogram(h, x0=b['low'][0], binsize=binsize, width=2)
        ph.label = 'All BOSS'

        #allhist.append(ph)
        #plt.add(ph)

        fill=biggles.FillBelow(x, h, color='black')
        fill.label = 'ALL BOSS'
        plt.add(fill)
        allhist.append(fill)

        """
        allmatches = self.read_matched('all')
        btot=eu.stat.Binner(allmatches['psf_fwhm_r'])
        btot.dohist(binsize=binsize, min=xmin, max=xmax)
        btot.calc_stats()
        htot = b['hist']/float(b['hist'].sum())
        """

        # first get the normalizations
        counts=0
        for i in xrange(ntype):
            type = types[i]
            t= self.read_matched(type)
            counts += t.size
            del t

        for i in xrange(ntype):
            type = types[i]

            t= self.read_matched(type)

            if 'psf_fwhm' in t.dtype.names:
                b=eu.stat.Binner(t['psf_fwhm'][:,2])
            else:
                b=eu.stat.Binner(t['psf_fwhm_r'])

            b.dohist(binsize=binsize, min=xmin, max=xmax)
            b.calc_stats()
            htot += b['hist']

            h = float(t.size)/counts*b['hist']/float(b['hist'].sum())
            #fracheight = float(t.size)/counts
            #h = fracheight*b['hist']/float(b['hist'].max())


            if name == 'all':
                if type == 'primus':
                    width=2
                else:
                    width=1
            else:
                width=2
            ph = biggles.Histogram(h, x0=b['low'][0], binsize=binsize, color=colors[i], width=width)
            ph.label = type

            allhist.append(ph)
            plt.add(ph)

        #htot = 0.8*htot/float(htot.max())
        htot = htot/float(htot.sum())
        ph = biggles.Histogram(htot, x0=b['low'][0], binsize=binsize, width=3, color=colors[-1])
        ph.label = 'Total'

        allhist.append(ph)
        plt.add(ph)

        """
        ph = biggles.Histogram(htot, x0=btot['low'][0], binsize=binsize, width=3, color='magenta')
        ph.label = 'All Matches'
        allhist.append(ph)
        plt.add(ph)
        """




        key=biggles.PlotKey(0.1, 0.95, allhist)

        plt.add(key)

        #plt.xrange = [0.4, 1.8]
        #if yrange is None:
        #    yrange = [0.0, 0.45]
        #yrange = [0.0, 1]
        #plt.yrange = yrange
        plt.xlabel = 'seeing FWHM [arcsec]'
        plt.aspect_ratio=0.7
        plt.show()

        epsfile=self.seeing_plotfile(name)
        print("Writing eps file:",epsfile)
        plt.write_eps(epsfile)
        converter.convert(epsfile, dpi=120, verbose=True)


