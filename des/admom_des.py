from sys import stderr
import numpy
from numpy import sqrt, zeros, ones
import esutil as eu
from esutil.numpy_util import where1
from esutil.misc import wlog
import fitsio
import admom
from . import collate

def process_image_cat(imfile, catfile, get_meta=False):
    """
    Process the input image and cat files.

    If get_meta=True, a tuple of (array, meta) is returned.  meta is a
    dictionary with 'shiftmax' and 'nsub'
    """
    im, ivar_im, cat = read_image_cat(imfile, catfile)

    col=cat['x_image']-1
    row=cat['y_image']-1
    sky=cat['background']

    colint=col.astype('i4')
    rowint=row.astype('i4')
    ivar = ivar_im[rowint,colint] 
    sigsky = zeros(cat.size) 
    w=where1(ivar > 0)
    wbad=where1(ivar <= 0)

    if w.size > 0:
        sigsky[w] = 1.0/sqrt(ivar[w])

        print >>stderr,'running admom'
        tmp=admom.admom(im, row[w], col[w], sky=sky[w], sigsky=sigsky[w])
        print >>stderr,'copying to array output'
        a,h=as_array(cat, w, wbad, tmp)
        if get_meta:
            out=(a,h)
        else:
            out=a
    else:
        out=None

    return out

def read_image_cat(imfile, catfile):
    im_ext = 1
    ivar_ext = 3
    cat_ext = 2

    if imfile[0:5] == 'hdfs':
        with eu.hdfs.HDFSFile(imfile,verbose=True) as fobj:
            fobj.stage()
            f = fitsio.FITS(fobj.localfile)
            im = f[im_ext][:,:]
            ivar_im = f[ivar_ext][:,:]
    else:
        with fitsio.FITS(imfile) as ff:
            im = ff[im_ext][:,:]
            ivar_im = ff[ivar_ext][:,:]

    cat = eu.io.read(catfile, type='fits', verbose=True, 
                     lower=True, ext=cat_ext)

    return im, ivar_im, cat

def as_array(cat, w, wbad, odict):
    dt=[]
    for n in ['row','col','Irr','Irc','Icc','e1','e2',
              'rho4','a4','s2','uncer','s2n','numiter',
              'wrow','wcol','whyflag','whystr',
              'sky','sigsky','guess']:
        if n in odict:
            dt += [(n,odict[n].dtype.str)]
    for n in cat.dtype.names:
        # only copy scalar fields
        if len(cat[n].shape) == 1:
            dt += [(n,cat[n].dtype.str)]

    arr=zeros(cat.size, dtype=dt) 

    eu.numpy_util.copy_fields(cat, arr)

    for n in odict:
        if n in arr.dtype.names:
            arr[n][w] = odict[n]

    if wbad.size > 0:
        arr['whyflag'][wbad] = 2**31
        arr['whystr'][wbad] = 'badivar'
        arr['Irr'][wbad] = -9999.
        arr['Irc'][wbad] = -9999.
        arr['Icc'][wbad] = -9999.
        arr['e1'][wbad] = -9999.
        arr['e2'][wbad] = -9999.
        arr['rho4'][wbad] = -9999.
        arr['uncer'][wbad] = 9999.
        arr['wrow'][wbad] = -9999.
        arr['wcol'][wbad] = -9999.
        arr['a4'][wbad] = -9999.
        arr['s2n'][wbad] = -9999.

    h={'shiftmax':odict['shiftmax'],'nsub':odict['nsub']}
    return arr, h


def compare_s2n(amrun, serun, matches=None):
    import biggles
    amc = collate.open_columns(amrun)
    sec = collate.open_columns(serun)

    if matches is None:
        wlog('reading am rid')
        arid = amc['rid'][:]
        wlog('  read:  ',arid.size)
        wlog('  unique:',numpy.unique(arid).size)

        wlog('reading se rid')
        srid = sec['rid'][:]
        wlog('  read:  ',srid.size)
        wlog('  unique:',numpy.unique(srid).size)

        wlog('matching rids')
        ma, ms = eu.numpy_util.match(arid, srid)

        if ma.size != arid.size:
            raise RuntimeError("matched only %d/%d" % (ma.size, arid.size))

    else:
        ma = matches['ma']
        ms = matches['ms']

    wlog("Reading am s2n, whyflag, Irr, Icc")
    am_s2n = amc['s2n'][:]
    am_irr = amc['Irr'][:]
    am_icc = amc['Icc'][:]
    #whyflag = amc['whyflag'][:]

    wlog("Reading se s2n,shear_flags")
    se_s2n = sec['shear_s2n'][:]
    seflags = sec['shear_flags'][:]

    #w=where1( (whyflag[ma] == 0) & (seflags[ms] == 0) )
    """
    w=where1(  (am_s2n[ma] > 0) 
             & (am_irr[ma] > 0) 
             & (am_icc[ma] > 0) 
             & (seflags[ms] == 0) )
    """
    w=where1(  (am_s2n[ma] > 0) & (seflags[ms] == 0) )

    msk = ms[w]
    mak = ma[w]

    wlog("min 'good' am s2n: ",am_s2n[mak].min())

    nperbin=100000

    wlog("binning nperbin:",nperbin)

    bs = eu.stat.Binner(se_s2n[msk], am_s2n[mak])
    bs.dohist(nperbin=nperbin, min=20.0, max=1000.0)
    bs.calc_stats()

    wlog("plotting")
    plt=eu.plotting.bscatter(bs['xmean'],bs['ymean'],yerr=bs['yerr'],
                             show=False,
                             xlabel=r'$s2n_{SH}$',ylabel=r'$s2n_{AM}$')


    coeff = numpy.polyfit(bs['xmean'], bs['ymean'], 1)
    poly=numpy.poly1d(coeff)

    flabt='m: %0.2f b: %0.3f' % (coeff[0],coeff[1])
    flab=biggles.PlotLabel(0.1,0.9,flabt,halign='left')
    plt.add(flab)


    ps = biggles.Curve(bs['xmean'], poly(bs['xmean']), color='blue')
    plt.add(ps)
    plt.show()

    plt.write_eps('/direct/astro+u/esheldon/tmp/compare-s2n.eps')

    wlog("binning sigma")
    sig=sqrt((am_irr[mak] + am_icc[mak])/2.)
    bs2 = eu.stat.Binner(sig, am_s2n[mak]/se_s2n[msk])
    #bs2.dohist(nperbin=nperbin)
    bs2.dohist(binsize=0.1, min=1, max=20)
    bs2.calc_stats()

    wlog("plotting")
    wp=where1(bs2['xmean'] > 0)
    plt2=eu.plotting.bscatter(bs2['xmean'][wp],
                              bs2['ymean'][wp],
                              yerr=bs2['yerr'][wp],
                              show=False,
                              xlabel=r'$\sigma_{AM} [pixels]$', 
                              ylabel=r'$s2n_{AM}/s2n_{SH}$')

    flatv = ones(bs2['xmean'][wp].size)*coeff[0]
    flat = biggles.Curve(bs2['xmean'][wp], flatv, color='blue')
    plt2.add(flat)
    plt2.show()
    plt2.write_eps('/direct/astro+u/esheldon/tmp/compare-s2n-vs-sigma.eps')

    return {'ma':ma, 'ms':ms}
