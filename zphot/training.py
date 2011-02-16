"""
See docs in weighting.py
"""
import os
from sys import stdout,stderr
import pprint
import copy

import numpy
from numpy import where

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

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
        #self.types = ['deep2','other','sdssobjids','vvds','zcosmos']
        self.types = ['primus.zerod.10oct29.zconf4',
                      '2slaq',
                      'cfrs',
                      'cnoc2.cut',
                      'deep2.form.fix',
                      'sdssobjids',
                      'tkrs-fix',
                      'vvds',
                      'zcosmos']
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
        dt.append( ('modelmag_dered','f4',5) )
        dt.append( ('modelmag_dered_err','f4',5) )
        dt.append( ('cmodelmag_dered','f4',5) )
        dt.append( ('cmodelmag_dered_err','f4',5) )

        if 'primus' in type:
            dt.append( ('z','f4') )
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

        if self.no_photo_cuts:
            procrun=self.photo_conf['procrun']
            import sdssgal
            cols = sdssgal.open_columns(procrun)

            stdout.write("\nReading ALL photoid,mag,cmag\n")
            photoid = cols['photoid'][:]
            mag = cols['modelmag_dered'][:]
            magerr = cols['modelmag_dered_err'][:]
            cmag = cols['cmodelmag_dered'][:]
            cmagerr = cols['cmodelmag_dered_err'][:]

            stdout.write('Reading ALL ra/dec/htmid\n')
            ra = cols['ra'][:]
            dec = cols['dec'][:]
            htmid = cols['htmid10'][:]
        else:
            photo_sample = self.conf["photo_sample"]
            zcs = zphot.select.ColumnSelector(photo_sample)

            zcs.select()

            indices = zcs.keep_indices.copy()

            stdout.write("\nReading photoid,mag,cmag\n")
            photoid = zcs.cols['photoid'][indices]
            mag = zcs.cols['modelmag_dered'][indices]
            magerr = zcs.cols['modelmag_dered_err'][indices]
            cmag = zcs.cols['cmodelmag_dered'][indices]
            cmagerr = zcs.cols['cmodelmag_dered_err'][indices]

            rmin = zcs.conf['cmodel_rmin']
            rmax = zcs.conf['cmodel_rmax']

            if (cmag[:,2].min() < rmin) or (cmag[:,2].max() > rmax):
                raise ValueError("Expected max rmag within [%s,%s]\n" % (rmin,rmax))

            # get all ra/dec for matching
            stdout.write('Reading ra/dec/htmid\n')
            ra = zcs.cols['ra'][indices]
            dec = zcs.cols['dec'][indices]
            htmid = zcs.cols['htmid10'][indices]

        h = eu.htm.HTM(10)

        minid = htmid.min()
        maxid = htmid.max()
        stdout.write("Getting htmrev\n")
        hist,htmrev = eu.stat.histogram(htmid-minid, rev=True)

        # match within 2 arcsec
        matchrad = 2.0/3600.0
        
        for type in self.types:
            stdout.write('\n' + '-'*70+'\n')
            spec = self.read_original(type)

            w = self.make_specific_cuts(spec, type)
            spec = spec[w]

            stdout.write("Matching\n")
            mspec,mphot,d12 = h.match(spec['ra'],spec['dec'], ra,dec,
                                      matchrad,
                                      htmid2 = htmid, 
                                      minid=minid,
                                      maxid=maxid,
                                      htmrev2=htmrev)

            mphot_u = numpy.unique1d(mphot)
            print '%s/%s are unique' % (mphot_u.size,mphot.size)
            stdout.write("Matched: %s/%s\n" % (mspec.size,spec.size))
            dt = self.matched_dtype(spec.dtype.descr, type)
            output = numpy.zeros(mspec.size,dtype=dt)
            eu.numpy_util.copy_fields(spec[mspec],output)

            if 'primus' in type:
                output['z'] = output['zprimus']
            
            stdout.write("copying photoid\n")
            output['photoid'] = photoid[mphot]

            stdout.write("copying mags\n")

            # note we don't want to read the matches directly from the columns
            # database this way because there are duplicates in mphot.  We must
            # keep in memory to do duplicate subscripting

            output['modelmag_dered']      = mag[mphot]
            output['modelmag_dered_err']  = magerr[mphot]
            output['cmodelmag_dered']     = cmag[mphot]
            output['cmodelmag_dered_err'] = cmagerr[mphot]

            # sanity check
            if not self.no_photo_cuts:
                if ((output['cmodelmag_dered'][:,2].max() > rmax) or
                    (output['cmodelmag_dered'][:,2].min() < rmin) ):
                    raise ValueError("Expected output rmag in [%s,%s]\n" % (rmin,rmax))

            # during matching we might end up with some too-faint mags for
            # the sdss
            if type == 'sdssobjids':
                stdout.write("    limiting SDSS training set to r < %s\n" \
                             % self.conf['sdss_rmax'])
                w,=where( output['cmodelmag_dered'][:,2] < self.conf['sdss_rmax'])
                stdout.write("        keeping %s/%s\n" % (w.size,mphot.size))
                output = output[w]
            

            # this may be sample dependant.
            zmin = self.conf['zmin']
            zmax = self.conf['zmax']
            stdout.write("Selecting spec to be in redshift "
                         "range: [%s,%s]\n" % (zmin,zmax))
            z = output['z']
            w,=where((z > zmin) & (z < zmax))
            stdout.write("Keeping %s/%s\n" % (w.size,mphot.size))
            output = output[w]



            outname=self.fname_matched(type)
            stdout.write("Writing to match file: %s\n" % outname)
            hdr = copy.deepcopy(self.conf)

            origfile=self.fname_original(type)
            hdr['original_file'] = os.path.basename(origfile)
            eu.io.write(outname, output, delim=' ', header=hdr)

    def make_specific_cuts(self, data, type):
        if type == 'deep2.form.fix':
            w,=where(data['zqual'] >= 3)
            stdout.write("deep2: keeping %s/%s\n" % (w.size, data.size))
            return w
        elif type == 'zcosmos':
            cc = data['cc']
            w,=where((cc==3.4) | (cc==3.5) | (cc==4.4) | (cc==4.5) | (cc==9.5))
            stdout.write("zcosmos: keeping %s/%s\n" % (w.size, data.size))
            return w
        elif type == 'vvds':
            zflag = data['zflag']
            w,=where((zflag==3) | (zflag==4))
            stdout.write("vvds: keeping %s/%s\n" % (w.size, data.size))
            return w
        else:
            return numpy.arange(data.size,dtype='i4')



    def read_original(self, type):
        fname = self.fname_original(type)
        stdout.write('Reading original spec file: %s\n' % fname)
        return eu.io.read(fname)

    def read_matched(self, type):
        fname = self.fname_matched(type)
        stdout.write('Reading matched spec file: %s\n' % fname)
        return eu.io.read(fname)


    def fname_original(self, type):
        ext='.rec'
        dir = path_join(self.basedir, 'original-2010-09-19')
        return os.path.join(dir, type+ext)

    def fname_matched(self, type):
        dir = self.dir_matched
        fname = type+'-match-'+self.train_sample+'.rec'
        fname = os.path.join(dir, fname)
        return fname

    def plotdir(self):
        dir = self.dir_matched
        return path_join(dir, 'plots')

    def plotfile(self, type):
        dir = self.plotdir()
        fname = self.fname_matched(type)
        fname = os.path.basename(fname)
        fname = fname.replace('.rec','.eps')
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

        stdout.write('\n')
        dir=self.plotdir()
        if not os.path.exists(dir):
            os.makedirs(dir)

        psfile = self.plotfile(type)

        orig = self.read_original(type)
        mat = self.read_matched(type)

        plt=FramedPlot() 

        symsize = 2
        if type == 'sdssobjids' or type == 'other':
            symsize = 0.25

        plt.add( Points(orig['ra'],orig['dec'],type='dot', size=symsize) )

        plt.add( Points(mat['ra'],mat['dec'],type='dot',
                        color='red', size=symsize) )
        plt.xlabel = 'RA'
        plt.ylabel = 'DEC'
        stdout.write("Writing eps file: %s\n" % psfile)
        plt.write_eps(psfile)
        converter.convert(psfile, dpi=120, verbose=True)


