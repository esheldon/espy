import os
from sys import stdout,stderr
import copy

import zphot
import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1
import sdssgal

import sdsspy

import numpy


bounds={}

bounds['zcosmos']   = {'ra': [149.44, 150.77], 'dec': [1.60,  2.85]}
bounds['primus1']   = {'ra': [36.5,     38.5], 'dec': [0.35,   0.9]}
bounds['primus2']   = {'ra': [309.5,   316.5], 'dec': [-0.9,   0.3]}
bounds['primus3']   = {'ra': [351.3,   353.6], 'dec': [-0.1,  0.44]}
bounds['deep2-egs'] = {'ra': [213.4,   215.6], 'dec': [51.9, 53.55]}
bounds['vvds']      = {'ra': [333.4,   335.5], 'dec': [-0.6,   1.3]}

def plot_bounds(doboss=False, **keys):
    """
    If procrun is set, a random subset of galaxies in the BOSS is 
    plotted in black, and the ones i the BOSS footprint are in
    grey
    """
    import biggles
    plt = biggles.FramedPlot()

    curves = []
    if doboss:
        f='~/masks/stomp-sdss/boss_survey.par'
        g = sdsspy.yanny.readone(f)
        i=0
        for b in g:
            # create curves in ra/dec space
            ra,dec = eu.plotting.transform_box(b['clambdaMin'],
                                                  b['clambdaMax'],
                                                  b['cetaMin'],
                                                  b['cetaMax'],
                                                  'sdss','eq')
            ra = eu.coords.shiftlon(ra,90)
            curve = biggles.Curve(ra,dec,color='grey')
            if i == 0:
                curve.label = 'BOSS'
                curves.append(curve)
            
            plt.add(curve)
            i+=1


    colors=['red','cyan','cadetblue','blue','magenta','orange']
    i=0
    for t in sorted(bounds):
        b=bounds[t]
        ra0 = eu.coords.shiftlon(b['ra'][0], 90)[0]
        ra1 = eu.coords.shiftlon(b['ra'][1], 90)[0]
        box = eu.plotting.bbox(ra0,ra1,b['dec'][0],b['dec'][1],
                               color=colors[i], linewidth=3)
        box.label = t
        plt.add(box)
        curves.append(box)
        i+=1

    key = biggles.PlotKey(0.7,0.95, curves)
    plt.add(key)
    plt.xlabel = 'RA - 90'
    #plt.xlabel = 'RA'
    plt.ylabel = 'DEC'

    epsfile = keys.get('epsfile',None)
    if epsfile is None:
        plt.show()
    else:
        epsfile = os.path.expandvars(epsfile)
        epsfile = os.path.expanduser(epsfile)
        plt.write_eps(epsfile)
        import converter
        converter.convert(epsfile,dpi=100,verbose=True)

class RachelQA:
    def __init__(self, pzrun, nchunk=50):
        """
        Take objects from the given pzrun and select:
            1) ones that are "recoverable"
            2) lie within the bounds of one of our test
            regions; these overlap spectroscopic regions.
        """

        self.pzrun=pzrun
        self.nchunk=nchunk
        self.conf = zphot.cascade_config(pzrun)
        self.wo = zphot.weighting.WeightedOutputs()

        # load columns
        self.procrun = self.conf['photo']['procrun']
        self.cols = sdssgal.open_columns(self.procrun)


    def process(self):
        """
        Run over all chunks, matching each against the
        bounds of each survey.  We will append to a file
        corresponding to each survey.
        """

        stdout.write("Reading all ra,dec from '%s'\n" % self.procrun)
        alldata = self.cols.read_columns(['ra','dec'])
        ra = self.cols['ra'][:]
        dec = self.cols['dec'][:]
        stdout.write("Read %s objects\n" % ra.size)


        # now get which ones are within the various bounds
        # this will remove most of the objects
        w = self.check_bounds(ra,dec)
        del ra
        del dec

        alldata = self.cols.read_columns(rows=w)

        num = self.wo.read_num(self.conf['weights']['wrun'], 1)

        r = self.open_output()

        num_beg = 0

        ntot = 0
        for chunk in xrange(self.nchunk):
        #for chunk in [0]:
            stdout.write('-'*70 + '\n')
            pzdata = zphot.weighting.read_pofz_byrun(self.pzrun,chunk)

            num_end = num_beg+pzdata.size
            tid = num['photoid'][num_beg:num_end]
            wbad = where1( pzdata['photoid'] != tid)
            if wbad.size > 0:
                eu.misc.colprint(wbad,pzdata['photoid'][wbad], tid[wbad])
                raise ValueError("Found non-matching photoid %d/%d" % \
                                 (wbad.size,pzdata.size))

            # grab subset with num > 0
            stdout.write("Getting num > 0\n")
            wnumgood = where1( num['num'][num_beg:num_end] > 0)
            stdout.write("    Found %s/%s\n" % (wnumgood.size, pzdata.size))

            # don't forget to increment
            num_beg += pzdata.size

            pzdata = pzdata[wnumgood]


            stdout.write("matching photoid\n")
            mall, mpz = eu.numpy_util.match(alldata['photoid'], pzdata['photoid'])

            if mall.size == 0:
                stdout.write("    NONE FOUND for chunk: %s\n" % chunk)
            else:
                ntot += mall.size
                stdout.write("    %s found for chunk: %s\n" % (mall.size,chunk))

                output = self.make_output(alldata, pzdata, mall, mpz)
                del pzdata
                del mall
                del mpz

                stdout.write("    Writing\n")
                r.write(output)
                del output
        stdout.write("Found total of %s matches\n" % ntot)
        r.close()

    def make_output(self, alldata, pzdata, mall, mpz):
        dt = self.out_dtype(alldata.dtype.descr)
        output = numpy.zeros(mall.size,dtype=dt)

        for n in alldata.dtype.names:
            output[n] = alldata[n][mall]

        #output['photoid'] = alldata['photoid'][mall]
        #output['ra'] = alldata['ra'][mall]
        #output['dec'] = alldata['dec'][mall]
        output['pofz'] = pzdata['pofz'][mpz]

        return output

    def out_dtype(self, old_dtype):
        nz = self.conf['pofz']['nz']
        #dt = [('photoid','i8'),('ra','f8'),('dec','f8'),
        #      ('pofz','f4',nz)]
        new_dtype = copy.deepcopy(old_dtype)
        new_dtype.append(('pofz','f4',nz))
        return new_dtype

    def outdir(self):
        outdir = self.wo.pofz_dir(self.pzrun)
        outdir = path_join(outdir, 'lenstest-match')
        return outdir

    def output_filename(self):
        fname = 'lenstest-match-'+self.pzrun+'.rec'
        fname = path_join(self.outdir(), fname)
        return fname

    def open_output(self):
        dir = self.outdir()
        if not os.path.exists(dir):
            os.makedirs(dir)
        fname = self.output_filename()
        stdout.write("Opening output file: '%s'\n" % fname)
        #r=eu.sfile.Open(fname,'w')
        r=eu.sfile.Open(fname,'w',delim=' ')
        return r

    def check_bounds(self, ra, dec):
        """
        Select all points within the bounds of the sample specified
        in "type"
        """

        logic = numpy.zeros(ra.size, dtype=numpy.bool)

        for type in sorted(bounds):
            stdout.write("Checking '%s'\n" % type)
            b = bounds[type]

            temp_logic = ( (ra >= b['ra'][0])   & (ra <= b['ra'][1]) 
                          & (dec >= b['dec'][0]) & (dec <= b['dec'][1]) )

            # or the logics together
            logic = logic | temp_logic
            w=where1(temp_logic)
            stdout.write("    Found %s\n" % w.size)

            del w

        wtot = where1(logic)
        stdout.write("Kept a total of %s\n" % wtot.size)
        return wtot
