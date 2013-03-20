import columns
import esutil as eu
from esutil.ostools import path_join
import biggles
import numpy
from numpy import where
import deswl

import sys
import os

if len(sys.argv) < 2:
    print 'compare-sizes.py run [show on screen]'
    sys.exit(45)

run=sys.argv[1]

cols=deswl.files.coldir_open(run)

collate_dir= deswl.files.wlse_collated_dir(run)
plotdir=path_join(collate_dir,'plots')



if not os.path.exists(plotdir):
    os.makedirs(plotdir)

band='i'
plotfile = os.path.join( plotdir, '%s-%s-sigma0-vs-shapelets-sigma.eps' % (run,band))


print 'reading data'
size_flags = cols['size_flags'][:]
shear_flags = cols['shear_flags'][:]
star_flag = cols['star_flag'][:]
ssigma=cols['shapelets_sigma'][:]
sigma0=cols['sigma0'][:]
psf_sigma=cols['interp_psf_sigma'][:]

print 'getting diff and sub-selecting good ones'
pdiff = numpy.sqrt( sigma0**2 - psf_sigma**2 )

w,=where( (size_flags == 0) & (shear_flags == 0) )

print 'histogramming'
ssh=eu.stat.histogram( ssigma[w], binsize=0.01, min=0, max=2)
s0h=eu.stat.histogram( sigma0[w], binsize=0.01, min=0, max=2)
diff_min = -0.1
pdiffh = eu.stat.histogram( pdiff[w], binsize=0.01, min=diff_min, max=2)


print 'making plot'
plt=biggles.FramedPlot()
sigma0_hist = biggles.Histogram(s0h, x0=0, binsize=0.01, color='blue')
sigma0_hist.label = 'sigma0'

ssigma_hist = biggles.Histogram(ssh, x0=0, binsize=0.01, color='red')
ssigma_hist.label='shapelets sigma'

pdiff_hist = biggles.Histogram( pdiffh, x0=diff_min, binsize=0.01, color='darkgreen')
pdiff_hist.label = 'sigma0 deconvolved'



key = biggles.PlotKey(0.55, 0.9, [sigma0_hist, ssigma_hist, pdiff_hist])

plt.add( sigma0_hist, ssigma_hist, pdiff_hist, key )

plt.xlabel = 'sigma'

if len(sys.argv) > 2:
    plt.show()

print 'writing to plot file:',plotfile
plt.write_eps( plotfile )

print 'converting to png'
eu.ostools.exec_process('converter -d 90 %s' % plotfile,verbose=True)


w2,=where( (size_flags == 0) & (shear_flags == 0) & (ssigma > 0.25) )

rat = sigma0[w2]/ssigma[w2]

wstar =where( (size_flags == 0) & (shear_flags == 0) & (star_flag == 1) )
star_rat = sigma0[wstar]/ssigma[wstar]

print 'mean ratio sigma0/ssigma:',rat.mean()
print 'mean ratio sigma0/ssigma for PSF stars:',star_rat.mean()
