from __future__ import print_function

import os
import numpy

def do_coadd_star_analysis(release, plot_sg=False, force=False):
    c=CoaddStarDensity(release.upper(), force=force)
    c.plot_num_hist()
    c.plot_radec()
    if plot_sg:
        c.plot_sg()

    return c

def get_dir():
    d=os.environ['DESDATA']
    d=os.path.join(d, 'users/esheldon/tile-star-density')
    return d

def get_sg_plotfile(release,tilename,band):
    d=get_dir()
    f='%s-%s-%s-sg.png' % (release,tilename,band)
    f=os.path.join(d,'sgplots',f)
    return f

def get_num_hist_plotfile(release):
    d=get_dir()
    f='%s-number-hist.png' % release
    f=os.path.join(d,f)
    return f

def get_redec_plotfile(release):
    d=get_dir()
    f='%s-radec.png' % release
    f=os.path.join(d,f)
    return f

def get_stardata_file(release):
    d=get_dir()
    f='%s-star-data.fits' % release
    f=os.path.join(d,f)
    return f



class CoaddStarDensity(object):
    """
    Measure the stellar density on all the coadds in the given
    release.  Make some plots.
    """
    def __init__(self, release, force=False, select_band='i'):
        self.select_band=select_band
        self.release=release

        self.nstar_cut=1700

        self._load_stardata(force=force)

        self._do_cut()

    def _do_cut(self):
        #w,=numpy.where(self.stardata['nstars'] > self.nstar_cut)
        w,=numpy.where(  ( (self.stardata['ra'] > 78)
                           & (self.stardata['dec'] < (-62)) )
                       | 
                         ( self.stardata['nstars'] > self.nstar_cut) )
        self.wbad=w

    def get_bad_tiles(self):
        bt=[t.strip() for t in self.stardata['tilename'][self.wbad]]
        return numpy.array(bt)

    def _load_stardata(self,force=False):
        import fitsio
        sdfile=get_stardata_file(self.release)
        if force or not os.path.exists(sdfile):
            self._load_runs()
            self._load_data()
            self._select_stars()
        else:
            self.stardata=fitsio.read(sdfile)


    def _load_runs(self):
        import desdb

        print("loading runs")
        self.coadd_runs = desdb.files.get_release_runs(self.release)
        self.run_data = self._get_run_data()
        #self.run_data = self.run_data[0:3]

    def _get_run_data(self):
        import desdb
        rlist= ','.join( ["'"+r+"'" for r  in self.coadd_runs] )

        query="""
        select
            distinct run,tilename
        from
            coadd
        where
            run in (%s)
        \n""" % rlist

        conn=desdb.Connection()
        res=conn.quick(query)
         
        return res

    def _load_data(self):
        import fitsio
        import desdb
        print("loading catalogs")

        df=desdb.DESFiles()

        alldata=[]

        cols=['flags','alphawin_j2000','deltawin_j2000','spread_model','mag_auto']

        nr=len(self.run_data)
        for i,rdata in enumerate(self.run_data):

            coadd_run=rdata['run']
            tilename=rdata['tilename']

            dd = {'coadd_run':coadd_run,
                  'tilename':tilename,
                  'band':self.select_band}

            f=df.url(type='coadd_cat',
                     coadd_run=coadd_run,
                     tilename=tilename,
                     band=self.select_band)

            print("    %d/%d  %s" % (i+1,nr,f))
            cat=fitsio.read(f,columns=cols,lower=True)

            dd['data'] = cat

            alldata.append(dd)
        
        self.alldata=alldata

    def _select_stars(self):
        import fitsio
        print("selecting stars")

        ntile=len(self.alldata)
        dt=[('tilename','S30'),
            ('ra','f8'),
            ('dec','f8'),
            ('nstars','i4')]

        stardata=numpy.zeros(ntile,dtype=dt)

        for i,dd in enumerate(self.alldata):
            cat=dd['data']

            w,=numpy.where(  (cat['flags']==0)
                           & (numpy.abs(cat['spread_model']) < 0.002)
                           & (cat['mag_auto'] < 19) & (cat['mag_auto'] > 16) )

            print("    ",dd['coadd_run'])
            print("        nstars:",w.size)
            dd['wstars'] = w

            stardata['tilename'][i] = dd['tilename']
            stardata['ra'][i] = numpy.median(cat['alphawin_j2000'][w])
            stardata['dec'][i] = numpy.median(cat['deltawin_j2000'][w])
            stardata['nstars'][i] = w.size

        self.stardata=stardata

        sdfile=get_stardata_file(self.release)
        print(sdfile)
        fitsio.write(sdfile,self.stardata,clobber=True)

    def plot_num_hist(self, show=False):
        import esutil as eu

        min_n=self.stardata['nstars'].min()
        max_n=self.stardata['nstars'].max()

        binsize=100
        plt0=eu.plotting.bhist(self.stardata['nstars'],
                               min=min_n,
                               max=max_n,
                               binsize=binsize,
                               show=False,
                               xlabel=r'$N_{star}$')

        w=self.wbad

        plt=eu.plotting.bhist(self.stardata['nstars'][w],
                              binsize=binsize,
                              min=min_n,
                              max=max_n,
                              show=show,
                              color='red',
                              plt=plt0)

        pngfile=get_num_hist_plotfile(self.release)
        print(pngfile)
        plt.write_img(800,800,pngfile)
        
    def plot_radec(self, show=False):
        import esutil as eu

        plt0=eu.plotting.bscatter(self.stardata['ra'],
                                  self.stardata['dec'],
                                  xlabel='RA',
                                  ylabel='DEC',
                                  show=False)

        w=self.wbad
        plt=eu.plotting.bscatter(self.stardata['ra'][w],
                                 self.stardata['dec'][w],
                                 show=show,
                                 color='red',
                                 plt=plt0)



        pngfile=get_redec_plotfile(self.release)
        print(pngfile)
        plt.write_img(800,800,pngfile)
 
    def plot_sg(self, show=False):
        import esutil as eu

        if not hasattr(self,'alldata'):
            self._load_stardata(force=True)

        for dd in self.alldata:
            coadd_run=dd['coadd_run']

            cat=dd['data']
            wstars=dd['wstars']

            w,=numpy.where(cat['flags']==0)
            plt0=eu.plotting.bscatter(cat['mag_auto'][w],
                                      cat['spread_model'][w],
                                      xrange=[15,22],
                                      yrange=[-0.01,0.01],
                                      show=False,
                                      type='filled circle',
                                      size=0.5)
            plt=eu.plotting.bscatter(cat['mag_auto'][wstars],
                                     cat['spread_model'][wstars],
                                     color='red',
                                     plt=plt0,
                                     type='filled circle',
                                     size=0.5,
                                     show=False)

            pngfile=get_sg_plotfile(self.release,dd['tilename'],dd['band'])
            print(pngfile)
            plt.write_img(800,800,pngfile)

            if show:
                plt.show()
                key=raw_input("hit a key (q to quit): ")
                if key=='q':
                    stop
