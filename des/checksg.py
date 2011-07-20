from __future__ import print_function
import os
import esutil as eu
from esutil.numpy_util import where1
import deswl

def plot_sg_all_exposures(run, ptype='png'):
    f=deswl.files.wlse_collated_path('wlse0011t', 'goodlist', ftype='json')
    flist = eu.io.read(f,verbose=True)
    for finfo in flist:
        expname = finfo['exposurename']
        print('-'*70)
        print("exposurename:",expname)
        plot_sg_by_exposure(run, expname, ptype=ptype)

def plot_sg_by_exposure(run, exposurename, ptype='png'):
    """
    This forces writing of png to avoid tons of windows
    """
    for ccd in xrange(1,62+1):
        plot_sg_by_ccd(run, exposurename, ccd, ptype=ptype)

def plot_sg_by_ccd(run, exposurename, ccd, ptype='screen'):
    import biggles
    #statfile = deswl.files.wlse_path(exposurename,ccd,'stat',serun=run)
    #stat = eu.io.read(statfile)

    #cat = eu.io.read(stat['cat'], ext=2, lower=True)
    #stars = eu.io.read(stat['output_files']['stars'])
    f = deswl.files.wlse_path(exposurename,ccd,'stars',serun=run)
    stars = eu.io.read(f)

    wstar = where1(stars['star_flag'] == 1)

    tab = biggles.Table(2,1)
    xrng = 8,17

    # spread model vs mag
    plt1 = biggles.FramedPlot()

    pall1 = biggles.Points(stars['mag'], stars['sg'], type='filled circle', size=0.25)
    pstar1 = biggles.Points(stars['mag'][wstar], stars['sg'][wstar], 
                            type='filled circle', size=0.85, color='red')

    plt1.add(pall1,pstar1)


    plt1.xrange = xrng
    plt1.yrange = -0.15,1.5
    plt1.ylabel = 'spread model'
    #plt1.x1.draw_ticklabels=0

    # sigma0 model vs mag
    plt2 = biggles.FramedPlot()

    pall2 = biggles.Points(stars['mag'], stars['sigma0'], type='filled circle', size=0.25)
    pstar2 = biggles.Points(stars['mag'][wstar], stars['sigma0'][wstar],
                            type='filled circle', size=0.85, color='red')

    plt2.add(pall2,pstar2)

    plt2.xrange = xrng
    plt2.yrange = 0.3,1.2
    plt2.ylabel = r'$\sigma_0$'
    plt2.xlabel = 'mag model'
    
    tab[0,0] = plt1
    tab[1,0] = plt2

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
            tab.show()
        else:
            subdir = ['checksg',exposurename]
            outf = deswl.files.wlse_test_path(run, subdir=subdir, extra='%02d' % ccd, ext=ext)
            d = os.path.dirname(outf)
            if not os.path.exists(d):
                print("making dir:",d)
                os.makedirs(d)
            print("Writing plot:",outf)
            if pt == 'eps':
                tab.write_eps(outf)
            else:
                tab.write_img(800,800,outf)


