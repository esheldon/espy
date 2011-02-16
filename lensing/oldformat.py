"""
This is primarily for dealing with the old outputs used in
the 2004 galaxy-mass correlation function paper
"""

import os
import sys
import numpy
import es_util as eu
import sfile

try:
    import pylab
    params = \
            {'axes.labelsize':16,
             'xtick.labelsize':14,
             'ytick.labelsize':14}
    pylab.rcParams.update(params)
except:
    sys.stderr.write('Could not import pylab\n')



olddir = os.path.expanduser('~/data/lensout/oldformat/comb')
lum_ext = '_N1.fit'
#lum_ext = '_N4.fit'

paper_stripes = '09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37'
#paper_stripes = '10_11_12_35_36_37'

def GetPaperDataLumFiles(band, redgal=False):

    red=''
    bin='num'
    if redgal:
        red='red'
        bin=''

    indir=os.path.join(olddir, 'sublum')
    indir=os.path.join(indir, band)

    f_front = os.path.join(indir, red+'lum')
    f_back = 'threebin'+bin+'_zgal_gal_stripe'+paper_stripes+ \
        '_gri_recorr_h_jack_comb_comoving'+lum_ext
    fnames = [f_front + num + f_back for num in ['1','2','3'] ]
    return fnames

def GetPaperDataLum(band, redgal=False):
    
    fnames = GetPaperDataLumFiles(band, redgal=redgal)
    dlist = [] 
    for f in fnames:
        sys.stdout.write('%s\n' % f)
        if os.path.exists(f):
            d,h = eu.mrdfits(f,1)
            dlist.append(d)
        else:
            raise ValueError('\tmissing file: %s' % f)
    data = eu.combine_arrlist(dlist)
    return data

def PlotPaperDataLum(band, redgal=False):
    d=GetPaperDataLum(band, redgal=redgal)

    #pylab.plot(d[0]['MEANR'],d[0]['SIGMA'],'.')
    pylab.clf()
    ax = pylab.subplot(111)

    ax.set_xscale('log')
    ax.set_yscale('log')
    for i in range(3):
        x=d[i]['MEANR_REBIN']
        y=d[i]['SIGMA_REBIN']
        yerr=d[i]['SIGMAERR_REBIN']
        ax.errorbar(x,y,yerr=yerr,fmt='.-')
    pylab.show()

def PaperDataLumToAscii(band, shortnames=False, redgal=False):
    import records
    d=GetPaperDataLum(band, redgal=redgal)
    fold = GetPaperDataLumFiles(band, redgal=redgal)

    red=''
    if redgal:
        red='red-'

    num=len(d)
    if shortnames:
        dirold=os.path.dirname(fold[0])
        fnames = [red+'deltasigma-'+band+'-'+str(n+1)+'.csv' for n in range(num)]
        fnames = [os.path.join(dirold,f) for f in fnames]
    else:
        fnames = [f.replace('.fit','.rec') for f in fold]

    dt=[('r','f4'),('delta_sigma','f4'),('delta_sigma_err','f4')] 
    for i in range(num):
        sys.stdout.write('%s\n' % fnames[i])
        nrad = len(d[i]['MEANR_REBIN'])
        out=numpy.zeros(nrad, dtype=dt)
        out['r'] = d[i]['MEANR_REBIN']
        out['delta_sigma'] = d[i]['SIGMA_REBIN']
        out['delta_sigma_err'] = d[i]['SIGMAERR_REBIN']
        
        if not shortnames:
            hdr={}
            hdr['file_origin'] = os.path.basename(fold[0])
        else:
            hdr=None
        sfile.write(out, fnames[i], header=hdr, delim=',') 
