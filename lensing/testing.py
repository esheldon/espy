import os
import numpy
from numpy import cos, sin, arccos, sqrt
import lensing
import esutil as eu
from esutil.coords import eq2xyz, gcirc, sphdist
from esutil.numpy_util import where1
from esutil.ostools import path_join

try:
    import biggles
    from biggles import FramedPlot,Points,Curve,PlotKey,PlotLabel
except:
    pass

def test_dsig(run, ntest, catdict=None, overwrite=False):

    outdir='/home/users/esheldon/www/tmp/test-shear'
    if catdict is None:
        catdict = load_test_data(run)

    conf = lensing.files.json_read(run)
    c = eu.cosmology.Cosmo(omega_m=conf['omega_m'])

    lcat = catdict['l']
    lout = catdict['lout']
    scat = catdict['s']

    xl=catdict['xl']
    yl=catdict['yl']
    zl=catdict['zl']
    xs=catdict['xs']
    ys=catdict['ys']
    zs=catdict['zs']

    m200_min = 3.e14
    wl = where1(lcat['m200'] > m200_min)

    rmin = 0.1
    #rmax = 1.0
    rmax = 0.5
    for ii in xrange(ntest):
        i = wl[ii]

        f = path_join(outdir,'lens-%04i.eps' % i)
        if not os.path.exists(f) or overwrite:
            print 'lens:',i
            print 'getting lens angular diameter distance'
            zlens   = lcat['z'][i]
            ralens  = lcat['ra'][i]
            declens = lcat['dec'][i]
            DL    = c.Da(0.0, zlens)[0]
            print '    zl: ',zlens
            print '    DL: ',DL

            maxangle = rmax/DL
            minangle = rmin/DL
            cos_maxangle = cos(maxangle)
            cos_minangle = cos(minangle)

            cosdis = xl[i]*xs + yl[i]*ys + zl[i]*zs

            w=where1( cosdis > cos_maxangle )

            print 'Number within %0.2f arcmin: %s/%s' \
                    % (maxangle*180/numpy.pi*60,w.size,scat.size)

            if w.size > 0:
                w2=where1( cosdis[w] < cos_minangle )
                if w2.size > 0:
                    w=w[w2]

                    # radians
                    r = arccos(cosdis[w])

                    # get total gamma instead bothering to project to tangent.
                    # Let's just see how it looks
                    #gamma = sqrt( s['g1'][w]**2 + s['g2'][w]**2 )
                    
                    # both radians
                    rr,theta = gcirc(ralens,
                                     declens,
                                     scat['ra'][w],
                                     scat['dec'][w],
                                     getangle=True)
                    rdiff = r-rr
                    print 'maxx diff rad:',numpy.abs(rdiff).max()
                    print theta[0:10]

                    r*=DL
                    rr*=DL
                    cos2theta = cos(2*theta)
                    sin2theta = sin(2*theta)
                    gammaT = -(scat['g1'][w]*cos2theta + scat['g2'][w]*sin2theta);

                    zsrc = scat['z'][w]

                    plt = FramedPlot()
                    p = Points(zsrc,gammaT*r,type='filled circle')
                    plt.add(p)
                    plt.xlabel = 'source z'
                    plt.ylabel = 'gamma*r (Mpc)'
                    plt.add(PlotLabel(0.2,0.9,'zlens: %0.2f' % zlens))
                    plt.add(PlotLabel(0.2,0.8,'m200: %0.2e' % lcat['m200'][i]))

                    print 'Writing to:',f
                    plt.write_eps(f)
                    #plt.show()

                    #tmp=raw_input('hit a key: ')
                    #if tmp == 'q':
                    #    return

def load_test_data(run):

    conf = lensing.files.json_read(run)

    s = lensing.files.scat_read(sample=conf['src_sample'])
    #l = lensing.files.lcat_read(sample=conf['lens_sample'])
    l = lensing.lcat.read_catalog(conf['lens_catalog'],conf['lens_version'])
    lout = lensing.files.lensout_read(run=run)

    print 'Getting lens x,y,z'
    xl,yl,zl = eq2xyz(l['ra'], l['dec'])
    print 'Getting source x,y,z'
    xs,ys,zs = eq2xyz(s['ra'], s['dec'])

    catdict = {}
    catdict['l'] = l
    catdict['lout'] = lout
    catdict['s'] = s

    catdict['xl'] = xl
    catdict['yl'] = yl
    catdict['zl'] = zl

    catdict['xs'] = xs
    catdict['ys'] = ys
    catdict['zs'] = zs

    return catdict



