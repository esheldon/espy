"""

Compare the princeton outputs with my regauss

"""

from __future__ import print_function
import os
import glob
import sdsspy
import numpy
import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join

import lensing

import biggles
from biggles import FramedPlot, PlotKey, Table, PlotLabel

def cache_all_matches(procrun,band):
    """
    Cache matching  for all runs
    """
    import pgnumpy
    pg = pgnumpy.PgNumpy()
    query = """
    select 
        distinct(run) 
    from 
        scat_princeton
    """
    print(query)
    runstr = pg.fetchall(query)
    runs = runstr['run']
    
    runs.sort()

    print("Found",len(runs)," runs")

    pc = Comparator(procrun,band)
    i=1
    for run in runs:
        print('-'*70)
        print("%d/%d" % (i,len(runs)))
        pc.load_matches(run, recache=True)
        if pc.match_data is not None:
            pc.write_matches(run)

        i+=1

class Comparator:
    def __init__(self,procrun, band):
        """

        For the input run, compare the new e1,e2 in pixel coords to the old
        princeton.  It *MUST* be on a run-by-run basis because of the orientations
        of runs.
        
        A problem with this is the princeton stuff was averaged r+i.

        read_matches() reads the run data from both princeton and mine and matches.
            If either data was missing for that run or no matches were found,
            .rgdata, .pgdata, etc. will be None

        write_matches(): write out the matches if found

        compare_run(): for input run, read data, do matches, and run various comparisions.
        compare_r(): Plot histograms for "r" for both data sets in the matching
            band.
        """

        import pgnumpy
        
        self.procrun=procrun
        self.band = band

        self.cols = lensing.regauss.open_columns(procrun)
        self.pg = pgnumpy.PgNumpy()
        self.current_run = None
        self.bandi = sdsspy.BANDS[band]

        self.reset()

    def reset(self):
        self.current_run=None
        self.pdata=None
        self.rgdata=None
        self.match_data = None

    def compare(self, run, recache=False):
        """
        run could be 'all'
        """
        self.load_matches(run, recache=recache)
        if self.match_data is None:
            return
        self.compare_r(run)
        self.compare_ellip(run)

    def compare_r(self, run, recache=False):
        self.load_matches(run, recache=recache)
        if self.match_data is None:
            return

        binsize=0.025
        min = 0.0
        max = 1.1

        tab = Table(2,1) 

        plt, pltdiff = self.compare_plots('r', min, max, binsize)
        tab[0,0] = plt
        tab[1,0] = pltdiff
        tab.show()


    def compare_ellip(self, run):
        self.load_matches(run, recache=recache)
        if self.match_data is None:
            return
        binsize=0.025
        min = -1.0
        max = 1.1

        tab = Table(2,2) 

        #
        # e1
        #

        column = 0
        for type in ['e1','e2']:

            plt, pltdiff = self.compare_plots(type, min, max, binsize)
            tab[0,column] = plt
            tab[1,column] = pltdiff

            column += 1

        tab.show()

    def compare_ellip_r(self):
        """
        Plot the median ellip as a function of r
        """
        pass

    def compare_plots(self, type, min, max, binsize):
        # full histograms
        pmed = numpy.median(self.match_data[type+'_p'])
        rgmed = numpy.median(self.match_data[type+'_rg'])

        plt,pph = eu.plotting.bhist(self.match_data[type+'_p'], 
                                    min=min,max=max,binsize=binsize, show=False,
                                    xlabel=type, 
                                    label='CH',getphist=True)

        plt,rgph = eu.plotting.bhist(self.match_data[type+'_rg'],
                                     min=min,max=max,binsize=binsize, show=False,
                                     plt=plt, 
                                     color='blue',
                                     label='rg '+self.band+'-band', getphist=True)

        key = PlotKey(0.1,0.9,[pph,rgph])
        namelab = PlotLabel(0.9,0.9,str(self.current_run), halign='left')

        pmedlab = PlotLabel(0.5,0.3,'p med: %0.2f' % pmed, halign='center')
        rgmedlab = PlotLabel(0.5,0.2,'rg med: %0.2f' % rgmed, halign='center')

        plt.add(key, namelab, pmedlab, rgmedlab)



        # histogram the differences
        diff = self.match_data[type+'_rg'] - self.match_data[type+'_p']
        print("mean",type,"diff:",diff.mean())
        print("median",type,"diff:",numpy.median(diff))

        pltdiff = eu.plotting.bhist(diff, binsize=binsize, show=False, 
                                    min=-2., max=2., 
                                    xlabel=r'%s$_{rg}$ - %s$_{CH}$' % (type,type))
        labmean = PlotLabel(0.1,0.9,'mean: %0.3f' % diff.mean(), halign='left')
        labmed = PlotLabel(0.1,0.8,'med:  %0.3f' % numpy.median(diff), halign='left')
        pltdiff.add(labmean,labmed)

        return plt, pltdiff


    def match_radec(self):
        matchrad_arcsec=1.0
        print("Matching by ra/dec within %0.1f arcsec" % matchrad_arcsec)
        h=eu.htm.HTM()
        mrg,mp,dis = h.match(self.rgdata['ra'],self.rgdata['dec'],
                             self.pdata['ra'],self.pdata['dec'],
                             1.0/3600.0)
        if mrg.size == 0:
            print("Found no ra/dec matches")
            return (None,None,None)

        print("Found",mrg.size,"matches")
        return mrg,mp,dis


    def load_matches(self, run, recache=False):
        """
        run could be 'all'

        Read data and reduce to the matches

        pdata and rgdata will be None if no data or matches are found
        """
        if run != self.current_run or recache:
            self.reset()
            self.current_run = run

            # try to read the cache
            if not recache:
                data = self.read_matches(run)
                if data is not None:
                    self.match_data = data
                    return

            print("Reading regauss data for run: ",run)

            if not self.cols['run'].has_index():
                self.cols['run'].create_index(db_verbose=1, tempdir='/dev/shm')
            wrun = (self.cols['run'] == run)
            if wrun.size == 0:
                self.reset()
                print("Run not found in regauss cols: %s" % run)
                return

            colnames = ['ra','dec','whyflag_rg','e1_rg','e2_rg','momerr_rg','r','r_rg','cmodelmag_dered']
            self.rgdata = self.cols.read_columns(colnames, rows=wrun)
            print("  Found",self.rgdata.size,"regauss objects")

            # minr to match princeton cut
            minr = 0.25
            w=where1(  (self.rgdata['r_rg'][:,self.bandi] > minr) 
                     & (self.rgdata['r_rg'][:,self.bandi] < 1.0)
                     & (self.rgdata['e1_rg'][:,self.bandi] > -4)
                     & (self.rgdata['e1_rg'][:,self.bandi] <  4)
                     & (self.rgdata['e2_rg'][:,self.bandi] > -4)
                     & (self.rgdata['e2_rg'][:,self.bandi] <  4)
                     & (self.rgdata['whyflag_rg'][:,self.bandi] == 0))
            if w.size == 0:
                self.reset()
                print("Found no object with r in [0.25,1], whyflag=0")
                return
            print("  Found",self.rgdata.size,"regauss objects with r in [0,1], whyflag=0")
            self.rgdata = self.rgdata[w]


            print("Reading princeton data for run: ",run)
            rname = 'r_'+self.band
            query="""
            select 
                ra,dec,e1pix as e1,e2pix as e2,e1e1err as momerr,{rname} as r
            from 
                scat_princeton
            where 
                run = {run}
                and {rname} > 0.0 
                and {rname} < 1.0
                and e1pix between -4.0 and 4.0
                and e2pix between -4.0 and 4.0
            """.format(run=run,rname=rname)
            print(query)
            self.pdata = self.pg.fetchall(query)

            if self.pdata.size == 0:
                self.reset()
                print("Run not found in princeton db: %s" % run)
                return
            print("  Found",self.pdata.size,"princeton objects")

            mrg,mp,dis = self.match_radec()
            if mrg is None:
                self.reset()
                return

            self.rgdata = self.rgdata[mrg]
            self.pdata = self.pdata[mp]

            self.copy_match_data()
            self.write_matches(self.current_run)

    def copy_match_data(self):
        if self.pdata is None:
            print("No matches are set")
            return
        st = self.match_struct(self.pdata.size)

        st['run']   = self.current_run
        st['band']  = self.band
        st['bandi'] = self.bandi
        st['e1_rg'] = self.rgdata['e1_rg'][:,self.bandi]
        st['e2_rg'] = self.rgdata['e2_rg'][:,self.bandi]
        st['r_rg']  = self.rgdata['r_rg'][:,self.bandi]

        st['e1_p']  = self.pdata['e1']
        st['e2_p']  = self.pdata['e2']
        st['r_p']   = self.pdata['r']

        st['rmag']  = self.rgdata['cmodelmag_dered'][:,2] 

        self.match_data = st

    def match_struct(self, n):
        dt=[('run',  'i2'),
            ('band', 'S1'),
            ('bandi','i1'),
            ('e1_rg','f4'),
            ('e2_rg','f4'),
            ('r_rg', 'f4'),
            ('e1_p', 'f4'),
            ('e2_p', 'f4'),
            ('r_p',  'f4'),
            ('rmag', 'f4')]
        return numpy.zeros(n, dtype=dt)

    def match_dir(self):
        dir='/home/users/esheldon/oh/lensinputs-v1/srcgal/princeton/match-my-regauss'
        return dir

    def match_file(self, run):
        """
        run might be 'all'
        """
        dir=self.match_dir()
        if not os.path.exists(dir):
            os.makedirs(dir)
        if run == 'all':
            f='match-regauss-%s-%s.fits' % (self.procrun,self.band)
        else:
            f='match-regauss-%s-%s-%06d.fits' % (self.procrun,self.band,run)

        f = path_join(dir,f)
        return f

    def read_matches(self, run):
        """
        run might be 'all'
        """
        f = self.match_file(run)
        if not os.path.exists(f):
            return None

        return eu.io.read(f,verbose=True)

    def write_matches(self, run):
        """
        run might be 'all'
        """
        if self.match_data is None:
            print("No matches are set")
            return

        f = self.match_file(run)
        eu.io.write(f, self.match_data, clobber=True, verbose=True)


    def combine_matches(self):
        """
        Combine all the caches into one big file.
        """
        dir = self.match_dir()
        pattern = path_join(dir,'match-regauss-*-*-*.fits')
        files = glob.glob(pattern)
        files.sort()
        print("Found",len(files),"files")

        data = eu.io.read(files, combine=True)
        outfile=self.match_file('all')
        eu.io.write(outfile,data,verbose=True)

