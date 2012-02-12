"""
    %prog [options] merun serun
"""
import os,sys
import recfile
import esutil as eu
from esutil.numpy_util import where1
import numpy
import deswl
import des
import converter

def make_html(pd,pfront):
    html="""
<html>
    <body bgcolor=white>

        <p>
        <img src="{pfront}-match-radiff-vs-shear-s2n.png">
        <img src="{pfront}-match-radiff-sdev-vs-shear-s2n.png">
        <p>
        <img src="{pfront}-match-decdiff-vs-shear-s2n.png">
        <img src="{pfront}-match-decdiff-sdev-vs-shear-s2n.png">

        <p>
        <img src="{pfront}-match-radiff-vs-mag-model.png">
        <img src="{pfront}-match-radiff-sdev-vs-mag-model.png">
        <p>
        <img src="{pfront}-match-decdiff-vs-mag-model.png">
        <img src="{pfront}-match-decdiff-sdev-vs-mag-model.png">
        <p>
        <img src="{pfront}-match-radiff.png">
        <img src="{pfront}-match-decdiff.png">

    </body>
</html>

    """.format(pfront=pfront)

    html_file=os.path.join(pd,pfront+'-poserr.html')
    print 'writing html file:',html_file
    with open(html_file,'w') as fobj:
        fobj.write(html)

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-s","--stars", action='store_true', help="Use PSF stars")

options, args = parser.parse_args(sys.argv[1:])


if len(args) < 2:
    parser.print_help()
    sys.exit(1)

merun=args[0]
serun=args[1]
d=deswl.files.collated_dir(merun)
pd=os.path.join(d,'poserr')
if not os.path.exists(pd):
    os.makedirs(pd) 

fextra=''
if options.stars:
    fextra='-stars'
    pfront='{serun}-stars-{merun}'.format(serun=serun,merun=merun)
else:
    pfront='{serun}-{merun}'.format(serun=serun,merun=merun)

f='{dir}/match-{serun}{fextra}-{merun}.dat'.format(dir=d,
                                                   merun=merun,serun=serun,
                                                   fextra=fextra)
f=os.path.expandvars(f)

print 'reading:',f
with recfile.Open(f,delim=' ',dtype=[('sid','i4'),('mid','i4')]) as fobj:
    matches = fobj.read()

print 'found:',matches.size,'matches'

mc = des.collate.open_columns(merun)
sc = des.collate.open_columns(serun)

mw = des.flags.select_good_me_bycol(mc)
if options.stars:
    star_flag=sc['star_flag'][:]
    sw=where1(star_flag != 0)
else:
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

make_html(pd,pfront)

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
        epsfile=pfront+'-match-radiff-sdev-vs-mag-model.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)

    if 'radiff_mean' in domag:
        print 'plotting mean radiff vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','radiff', 
                                   names=names,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile=pfront+'-match-radiff-vs-mag-model.eps'
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
        epsfile=pfront+'-match-decdiff-sdev-vs-mag-model.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'decdiff_mean' in domag:
        print 'plotting mean decdiff vs mag_model'
        plt = eu.plotting.bhist_vs(ddata, 'mag_model','decdiff', 
                                   names=names,
                                   max=25,
                                   nperbin=nperbin, xlog=xlog)
        epsfile=pfront+'-match-decdiff-vs-mag-model.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    del mag
    del magm


if ptypes['radiff_hist']:
    sdev=radiff.std()
    binsize = 0.1*sdev
    print 'radiff sdev:',sdev

    print 'doing radiff'
    plt=eu.plotting.bhist(radiff, binsize=binsize, xlabel=names['radiff'],ylog=True)
    epsfile=pfront+'-match-radiff.eps'
    epsfile=os.path.join(pd,epsfile)
    plt.write_eps(epsfile)
    converter.convert(epsfile, dpi=dpi, verbose=True)

if ptypes['decdiff_hist']:
    sdev=decdiff.std()
    binsize = 0.1*sdev
    print 'decdiff sdev:',sdev

    print 'doing decdiff'
    plt=eu.plotting.bhist(decdiff, binsize=binsize, xlabel=names['decdiff'],ylog=True)
    epsfile=pfront+'-match-decdiff.eps'
    epsfile=os.path.join(pd,epsfile)
    plt.write_eps(epsfile)
    converter.convert(epsfile, dpi=dpi, verbose=True)




if options.stars:
    sys.exit(0)


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
        epsfile=pfront+'-match-radiff-sdev-vs-shear-s2n.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'radiff_mean' in dos2n:
        print 'plotting mean radiff vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','radiff', 
                                   names=names,
                                   nperbin=nperbin, xlog=xlog)
        epsfile=pfront+'-match-radiff-vs-shear-s2n.eps'
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
        epsfile=pfront+'-match-decdiff-sdev-vs-shear-s2n.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


    if 'decdiff_mean' in dos2n:

        print 'plotting mean decdiff vs s2n'
        plt = eu.plotting.bhist_vs(ddata, 'shear_s2n','decdiff', 
                                   names=names,
                                   nperbin=nperbin, xlog=xlog)
        epsfile=pfront+'-match-decdiff-vs-shear-s2n.eps'
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)



    del s2n
    del s2nm


