import os
import sys
from numpy import where
import esutil as eu
import biggles
import pcolors
import glob
import converter
import cPickle

from . import files

def make_pointings(simname, frac=0.01):

    conf=files.read_config(simname)
    dir = files.get_pointings_dir(simname)
    purl = files.get_pointings_url(simname)
    eu.ostools.makedirs_fromfile(purl)

    print 'will write to:',purl

    pixscale_deg=conf['pixscale']/3600.

    t=files.read_original_list(conf['orig_vers'], fnum_list=conf['fields'])

    plt,ind=eu.plotting.plotrand(t['ra'], t['dec'],show=False,frac=frac,
                                 get_indices=True)

    ramin=t['ra'].min()
    ramax=t['ra'].max()
    decmin=t['dec'].min()
    decmax=t['dec'].max()

    #rra=t['ra']
    #rdec=t['dec']
    rra=t['ra'][ind]
    rdec=t['dec'][ind]

    nrow=4096
    ncol=4096
    rawidth=4096*pixscale_deg
    decwidth=4096*pixscale_deg

    ncenters_ra = 1 + int( (ramax-ramin)/rawidth )
    ncenters_dec = 1 + int( (decmax-decmin)/decwidth )

    ncolor=ncenters_ra*ncenters_dec
    colors=pcolors.rainbow(ncolor)
    
    hdict={}

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


            h={'crval1': ram,
               'crval2': decm,
               'crpix1': ncol/2.,
               'crpix2': nrow/2.,
               'cd1_1': 0.,
               'cd2_2': 0.,
               'ctype1': 'RA---TAN',
               'ctype2': 'DEC--TAN'}

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
                    pnum += 1
                    hdict[pnum] = h

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
    with open(purl,'w') as fobj:
        cPickle.dump(hdict, fobj)

    return


