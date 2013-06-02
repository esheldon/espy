"""
This algorithm is deeply flawed

I'm hitting memory issues too eve
"""
import os
import sys
import numpy
from numpy import where
import esutil as eu
import pcolors
import glob
import converter
import cPickle

from . import files

def get_struct():
    dt =   [('index','i4'),
            ('crval1','f8'),
            ('crval2','f8'),
            ('crpix1','f8'),
            ('crpix2', 'f8'),
            ('cd1_1','f8'),
            ('cd1_2','f8'),
            ('cd2_2','f8'),
            ('cd2_1','f8'),
            ('ctype1','S8'),
            ('ctype2','S8')]

    return numpy.zeros(1, dtype=dt)

def make_pointings(simname, frac=0.01):
    import biggles

    conf=files.read_config(simname)


    dir = files.get_pointings_dir(simname)
    purl = files.get_pointings_url(simname)
    eu.ostools.makedirs_fromfile(purl)

    print 'will write to:',purl

    pixscale_deg=conf['pixscale']/3600.

    # should use database instead
    t=files.read_original_list(conf['orig_vers'], fnum_list=conf['fields'])

    plt,ind=eu.plotting.plotrand(t['ra'], t['dec'],show=False,frac=frac,
                                 seed=3000,
                                 get_indices=True)

    ramin=t['ra'].min()
    ramax=t['ra'].max()
    decmin=t['dec'].min()
    decmax=t['dec'].max()

    #rra=t['ra']
    #rdec=t['dec']
    rra=t['ra'][ind]
    rdec=t['dec'][ind]

    nrow=conf['nrow']
    ncol=conf['ncol']

    # this is not correct!
    rawidth=ncol*pixscale_deg
    decwidth=nrow*pixscale_deg

    ncenters_ra = 1 + int( (ramax-ramin)/rawidth )
    ncenters_dec = 1 + int( (decmax-decmin)/decwidth )

    ncolor=ncenters_ra*ncenters_dec
    colors=pcolors.rainbow(ncolor)
    
    hlist=[]

    rastart=ramin + rawidth/2.
    decstart=decmin + decwidth/2.

    itot=0
    pnum=0
    for ira in xrange(ncenters_ra):
        ram=rastart + ira*rawidth
        print '%s/%s  %s' % (ira+1,ncenters_ra,ram)
        for idec in xrange(ncenters_dec):
            decm=decstart + idec*decwidth

            color=colors[itot % ncolor]
            itot += 1

            h=get_struct()

            h['crval1'] = ram
            h['crval2'] = decm
            h['crpix1'] = ncol/2.
            h['crpix2'] = nrow/2.
            h['cd1_1'] = 0.
            h['cd2_2'] = 0.
            h['ctype1'] = 'RA---TAN'
            h['ctype2'] = 'DEC--TAN'

            if conf['alt_cd']:
                h['cd1_2'] = pixscale_deg
                h['cd2_1'] = pixscale_deg
            else:
                h['cd1_2'] = pixscale_deg
                h['cd2_1'] = -pixscale_deg


            w,=where(  (rra  > (ram - rawidth/2 - 0.2)) 
                     & (rra  < (ram + rawidth/2 + 0.2)) 
                     & (rdec > (decm - decwidth/2 - 0.2)) 
                     & (rdec < (decm + decwidth/2 + 0.2)) )

            if w.size > 0:
                wcs=eu.wcsutil.WCS(h)
                col,row=wcs.sky2image(rra[w], rdec[w],find=False)
                w2,=where( (col > 0) 
                          &(col <= ncol)
                          &(row > 0)
                          &(row <= nrow) )
                if w2.size > 0:
                    rap,decp = wcs.image2sky([1,ncol],[1,nrow])
                    box=eu.plotting.bbox(rap.min(), rap.max(),
                                         decp.min(), decp.max(),
                                         width=0.75,color=color)
                    plt.add(box)

                    plt.add(biggles.PlotLabel(ram, decm,
                                              str(pnum),
                                              color=color,
                                              fontsize=0.5) )


                    h['index'] = pnum
                    hlist.append(h)

                    pnum += 1

    plt.xlabel='RA'
    plt.ylabel='DEC'
    plt.aspect_ratio=(decmax-decmin)/(ramax-ramin)

    epsfile=os.path.join(dir, 'pointings.eps')

    print epsfile
    plt.write_eps(epsfile)

    pngfile=epsfile.replace('.eps','.png')
    print pngfile
    plt.write_img(1200,1200,pngfile)

    print purl
    hdata=eu.numpy_util.combine_arrlist(hlist)
    eu.io.write(purl, hdata, clobber=True)

    return


