from __future__ import print_function
import os
from sys import stderr
import esutil as eu
from esutil.numpy_util import where1
import deswl
import numpy

class SizeMagPlotter(dict):
    def __init__(self, run, **keys):

        self['run'] = run
        for k,v in keys.iteritems():
            self[k] = v
        if 'ptype' not in self:
            self['ptype'] = 'png'
        if 'fs' not in self:
            self['fs'] = 'nfs'

        if self['ptype'] in ['png','eps']:
            self['fext']=self['ptype']
        else:
            self['fext']=None

        self._load_columns_if_exist()


    def plot_exposure(self, expname, combine=True):
        """
        Plot a size mag diagram for an exposure

        parameters
        ----------
        expname: string
            des single epoch exposure name
        combine: bool, optional
            Show a combined plot from all ccds, default False
        """

        import biggles
        import pcolors

        if not combine:
            # default to png here
            if self['ptype'] == 'screen':
                raise ValueError("don't use ptype == screen for plots by exposure "
                                 "unless you send combine=True also")
            for ccd in xrange(1,62+1):
                self.plot_ccd(expname, ccd)
        else:

            # make a combined plot
            colors = pcolors.rainbow(62, 'hex')
            tab=None
            keys={'doplot':False,'addlabel':False}
            for ccd in xrange(1,62+1):
                keys['tab'] = tab
                keys['star_color'] = colors[ccd-1] 
                tab=self.plot_ccd(expname, ccd, **keys)

            tab.title = '%s-%s' % (self['run'],expname)
            

            subdir = ['checksg','byexp']
            if self['ptype'] == 'screen':
                tab.show()
            else:
                outf = deswl.files.se_test_path(self['run'], 
                                                subdir=subdir, 
                                                extra=expname, fext=self['fext'])
                eu.ostools.makedirs_fromfile(outf)
                if self['ptype'] == 'eps':
                    print("writing plot file:",outf)
                    tab.write_eps(outf)
                else:
                    print("writing plot file:",outf)
                    tab.write_img(800,800,outf)



    def plot_ccd(self, expname, ccd, 
                 star_color='red', doplot=True, 
                 fwhm_range=None, tab=None, addlabel=True, **keys):
        import biggles

        if self.cols is None:
            f = deswl.files.se_url(expname,ccd,'stars',serun=self['run'],fs=self['fs'])
            data = eu.io.read(f)
            magname='mag'
        else:
            w=where1((self.expnames == expname) & (self.ccds == ccd))
            data=self.cols.read_columns(['star_flag','imag','sg','sigma0'],rows=w)
            magname='imag'

        wstar = where1(data['star_flag'] == 1)
        print('    expname:',expname,'ccd: %02i' % ccd,'Got',data.size,'objects,',wstar.size,'stars',file=stderr)

        xrng = 8,17

        # spread model vs mag
        if tab is None:
            tab = biggles.Table(2,1)
            plt1 = biggles.FramedPlot()
            plt2 = biggles.FramedPlot()
            tab[0,0] = plt1
            tab[1,0] = plt2

        if addlabel:
            lab = biggles.PlotLabel(0.1,0.9,'%s-%s-%0.2d' % (self['run'],expname,ccd), 
                                    halign='left')
            tab[0,0].add(lab)

        pall1 = biggles.Points(data[magname], data['sg'], type='filled circle', size=0.25)
        pstar1 = biggles.Points(data[magname][wstar], data['sg'][wstar], 
                                type='filled circle', size=0.85, color=star_color)

        tab[0,0].add(pall1,pstar1)


        tab[0,0].xrange = xrng
        tab[0,0].yrange = -0.01,0.03
        tab[0,0].ylabel = 'spread model'

        # fwhm model vs mag: sigma0 is in arcsec
        fwhm = 2.35*data['sigma0']
        pall2 = biggles.Points(data[magname], fwhm, type='filled circle', size=0.25)
        pstar2 = biggles.Points(data[magname][wstar], fwhm[wstar],
                                type='filled circle', size=0.85, color=star_color)

        tab[1,0].add(pall2,pstar2)

        tab[1,0].xrange = xrng
        if fwhm_range is None:
            tab[1,0].yrange = 0.6,2
        else:
            tab[1,0].yrange = fwhm_range
            
        #tab[1,0].ylabel = r'$\sigma_0$'
        tab[1,0].ylabel = 'FWHM'
        tab[1,0].xlabel = 'mag model'
        
        
        if doplot:
            if self['ptype'] == 'screen':
                tab.show()
            else:
                subdir = ['checksg',expname]
                outf = deswl.files.se_test_path(self['run'], 
                                                subdir=subdir, extra='%02d' % ccd, fext=self['fext'])
                eu.ostools.makedirs_fromfile(outf)
                print("Writing plot:",outf)
                if self['ptype'] == 'eps':
                    tab.write_eps(outf)
                else:
                    tab.write_img(800,800,outf)


        return tab



    def _load_columns_if_exist(self):
        self.cols=None
        return
        dir=deswl.files.coldir(self['run'])
        if not os.path.exists(dir):
            return
        self.cols=deswl.files.coldir_open(self['run'])


