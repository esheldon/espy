"""

Plot histograms of ra and dec diff between coadd and matched single epoch

"""

import os,sys
import recfile
import esutil as eu
from esutil.numpy_util import where1
import numpy
import des
import converter

def select_me(c):
    flags   = c['shear_flags'][:]
    flagsin = c['input_flags'][:]
    flagswt = c['flag_weights'][:]
    w = where1((flags == 0) & (flagsin == 0) & (flagswt==0))

    return w

def select_se(c):
    flags   = c['shear_flags'][:]
    w = where1(flags == 0)

    return w

if len(sys.argv) < 3:
    print 'usage: plot-match-posdiff.py merun serun'
    sys.exit(1)

merun=sys.argv[1]
serun=sys.argv[2]
d=deswl.files.collated_dir(merun)
pd=os.path.join(d,'poserr')
os.makedirs(pd) 

f='%(dir)s/match-%(serun)s-%(merun)s.dat' % {'dir':d,'merun':merun,'serun':serun}
f=os.path.expandvars(f)

print 'reading:',f
with recfile.Open(f,delim=' ',dtype=[('sid','i4'),('mid','i4')]) as fobj:
    matches = fobj.read()

print 'found:',matches.size,'matches'

mc = des.collate.open_columns(merun)
sc = des.collate.open_columns(serun)

mw = des.flags.select_good_me_bycol(mc)
sw = des.flags.select_good_se_bycol(sc)

print 'reading ra,dec from me'
mra = mc['ra'][mw]
mdec = mc['dec'][mw]

print 'reading ra,dec from se'
sra = sc['ra'][sw]
sdec = sc['dec'][sw]

dpi=90


ptypes = {'radiff_hist':True,
          'decdiff_hist':True,
          'mag':['radiff_mean','radiff_sdev','decdiff_mean','decdiff_sdev'],
          's2n':['radiff_mean','radiff_sdev','decdiff_mean','decdiff_sdev']}
#ptypes = {'radiff_hist':False,
#          'decdiff_hist':False,
#          'mag':False,
#          's2n':['radiff_sdev','decdiff_sdev']}

radiff = (mra[matches['mid']] - sra[matches['sid']])*3600
decdiff = (mdec[matches['mid']] - sdec[matches['sid']])*3600

names={'radiff':r'$\langle RA_{coadd}-RA_{se}\rangle$ [arcsec]',
       'decdiff':r'$\langle DEC_{coadd}-DEC_{se}\rangle$ [arcsec]'}
names_sdev = {'radiff':r'$\sigma(RA_{coadd}-RA_{se})$ [arcsec]',
              'decdiff':r'$\sigma(DEC_{coadd}-DEC_{se})$ [arcsec]'}

nperbin=100000

domag=ptypes['mag']
if domag:
    print 'reading mag_model'
    mag=mc['mag_model'][mw]

    magm = mag[matches['mid']]

    w,=numpy.where(magm > 0.1)

    ddata={'radiff':radiff[w], 'mag_model':magm[w]}

    xlog=False

    if 'radiff_sdev' in domag:
        print 'plotting std(radiff) vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','radiff', 
                                   stype='sdev',
                                   names=names_sdev,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_radiff_sdev_vs_mag_model.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)

    if 'radiff_mean' in domag:
        print 'plotting mean radiff vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','radiff', 
                                   names=names,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_radiff_vs_mag_model.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)

    ddata={'decdiff':decdiff[w], 'mag_model':magm[w]}

    if 'decdiff_sdev' in domag:
        print 'plotting std(decdiff) vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','decdiff', 
                                   stype='sdev',
                                   names=names_sdev,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_decdiff_sdev_vs_mag_model.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'decdiff_mean' in domag:
        print 'plotting mean decdiff vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','decdiff', 
                                   names=names,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_decdiff_vs_mag_model.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    del mag
    del magm


dos2n=ptypes['s2n']
if dos2n:
    print 'reading s2n'
    s2n=mc['shear_s2n'][mw]

    s2nm = s2n[matches['mid']]

    w,=numpy.where(s2nm > 0.1)
    ddata={'radiff':radiff[w], 'shear_s2n':s2nm[w]}

    xlog=True

    if 'radiff_sdev' in dos2n:
        print 'plotting std(radiff) vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','radiff', 
                                   stype='sdev',
                                   names=names_sdev,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_radiff_sdev_vs_shear_s2n.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'radiff_mean' in dos2n:
        print 'plotting mean radiff vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','radiff', 
                                   names=names,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_radiff_vs_shear_s2n.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)



    ddata={'decdiff':decdiff[w], 'shear_s2n':s2nm[w]}
    if 'decdiff_sdev' in dos2n:
        print 'plotting std(decdiff) vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','decdiff', 
                                   stype='sdev',
                                   names=names_sdev,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_decdiff_sdev_vs_shear_s2n.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'decdiff_mean' in dos2n:

        print 'plotting mean decdiff vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','decdiff', 
                                   names=names,
                                   nperbin=nperbin, xlog=xlog)
        epsfile='%s_%s_match_decdiff_vs_shear_s2n.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)



    del s2n
    del s2nm

if ptypes['radiff_hist']:
    sdev=radiff.std()
    binsize = 0.1*sdev
    print 'radiff sdev:',sdev

    if False:
        print 'doing radiff'
        plt=eu.plotting.bhist(radiff, binsize=binsize, xlabel=names['radiff'],ylog=True)
        epsfile='%s_%s_match_radiff.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt.write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)

if ptypes['decdiff_hist']:
    sdev=decdiff.std()
    binsize = 0.1*sdev
    print 'decdiff sdev:',sdev

    if False:
        print 'doing decdiff'
        plt=eu.plotting.bhist(decdiff, binsize=binsize, xlabel=names['decdiff'])
        epsfile='%s_%s_match_decdiff.eps' % (merun,serun)
        epsfile=os.path.join(pd,epsfile)
        plt.write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


