from __future__ import print_function
import os
import esutil as eu
from esutil.numpy_util import where1
import deswl

def plot_sg_all_exposures(run, **keys):
    if 'ptype' not in keys:
        keys['ptype'] = 'eps'

    f=deswl.files.wlse_collated_path('wlse0011t', 'goodlist', ftype='json')
    flist = eu.io.read(f,verbose=True)
    # get unique exposures
    d={}
    for f in flist:
        expname = f['exposurename']
        d[expname] = expname

    i=1
    n=len(d)
    for expname in sorted(d):
        print('-'*70)
        print("exposurename: %s  (%d/%d)" % (expname,i,n))
        plot_sg_by_exposure(run, expname, **keys)
        i+=1

def plot_sg_by_exposure(run, exposurename, **keys):
    """

    parameters
    ----------
    run: string
        serun
    exposurename: string
        des exposurename
    combine: bool, optional
        Show a combined plot from all ccds, default False

    And keywords to plot_sg_by_ccd, such as ptype=, fwhm_range=

    if not combine=True, this forces writing of png or eps to avoid tons of
    windows

    """

    import biggles

    combine = keys.get('combine',False)
    if not combine:
        # default to png here
        ptype = keys.get('ptype','png')
        if ptype == 'screen':
            raise ValueError("don't send ptype == screen for plots by exposure "
                             "unless you send combine=True also")
        for ccd in xrange(1,62+1):
            plot_sg_by_ccd(run, exposurename, ccd, **keys)
    else:

        # make a combined plot
        import pcolors
        colors = pcolors.rainbow(62, 'hex')
        tab=None
        for ccd in xrange(1,62+1):
            tab=plot_sg_by_ccd(run, exposurename, ccd, star_color=colors[ccd-1], 
                               tab=tab, show=False, addlabel=False)

        #lab = biggles.PlotLabel(0.05,0.9,'%s-%s' % (run,exposurename), halign='left',color='blue')
        #tab[0,0].add(lab)
        tab.title = '%s-%s' % (run,exposurename)
        
        subdir = 'checksg'
        ptype = keys.get('ptype','screen')
        if ptype == 'screen':
            tab.show()
        elif ptype == 'eps':
            outf = deswl.files.wlse_test_path(run, subdir=subdir, extra=exposurename, ext='eps')
            makedirs_fromfile(outf)
            print("writing plot file:",outf)
            tab.write_eps(outf)
        else:
            outf = deswl.files.wlse_test_path(run, subdir=subdir, extra=exposurename, ext='png')
            makedirs_fromfile(outf)
            print("writing plot file:",outf)
            tab.write_img(800,800,outf)


def plot_sg_by_ccd(run, exposurename, ccd, star_color='red', show=True, ptype='screen', fwhm_range=None, tab=None, addlabel=True):
    import biggles

    f = deswl.files.wlse_path(exposurename,ccd,'stars',serun=run)
    stars = eu.io.read(f)

    wstar = where1(stars['star_flag'] == 1)

    xrng = 8,17

    # spread model vs mag
    if tab is None:
        tab = biggles.Table(2,1)
        plt1 = biggles.FramedPlot()
        plt2 = biggles.FramedPlot()
        tab[0,0] = plt1
        tab[1,0] = plt2

    if addlabel:
        lab = biggles.PlotLabel(0.1,0.9,'%s-%s-%0.2d' % (run,exposurename,ccd), halign='left')
        tab[0,0].add(lab)

    pall1 = biggles.Points(stars['mag'], stars['sg'], type='filled circle', size=0.25)
    pstar1 = biggles.Points(stars['mag'][wstar], stars['sg'][wstar], 
                            type='filled circle', size=0.85, color=star_color)

    tab[0,0].add(pall1,pstar1)


    tab[0,0].xrange = xrng
    tab[0,0].yrange = -0.15,1.5
    tab[0,0].ylabel = 'spread model'

    # fwhm model vs mag: sigma0 is in arcsec
    fwhm = 2.35*stars['sigma0']
    pall2 = biggles.Points(stars['mag'], fwhm, type='filled circle', size=0.25)
    pstar2 = biggles.Points(stars['mag'][wstar], fwhm[wstar],
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
    

    if not isinstance(ptype,list):
        ptype = [ptype]

    for pt in ptype:
        if pt == 'eps':
            ext='eps'
        elif pt == 'png':
            ext='png'
        elif pt == 'screen':
            pass
        else:
            raise ValueError("each ptype should be 'screen','eps','png'")

        if pt == 'screen':
            if show:
                tab.show()
        else:
            subdir = ['checksg',exposurename]
            outf = deswl.files.wlse_test_path(run, subdir=subdir, extra='%02d' % ccd, ext=ext)
            makedirs_fromfile(outf)
            print("Writing plot:",outf)
            if pt == 'eps':
                tab.write_eps(outf)
            else:
                tab.write_img(800,800,outf)


    return tab

def makedirs_fromfile(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        print("making dir:",d)
        os.makedirs(d)

