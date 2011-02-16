import sys
from sys import stdout
import deswl
import columns
import biggles
import esutil
from esutil.ostools import path_join

if len(sys.argv) < 2:
    raise RuntimeError("usage: plot_imaghist serun")

serun=sys.argv[1]
coldir = deswl.files.wlse_coldir(serun)

c = columns.Columns(coldir)

w=c['psfstars']['psf_flags'] == 0

uid = c['psfstars']['uid'][w]

imag = c['imag'][uid]


plt=esutil.plotting.bhist(imag, binsize=0.01,xlabel='i-mag',show=False)

label=biggles.PlotLabel(.1, .9, "psf_flags==0",halign='left')
plt.add(label)


collated_dir = deswl.files.wlse_collated_dir(serun)
plotdir = path_join(collated_dir, 'plots')

plotfile='imaghist-%s.eps' % serun
plotfile=path_join(plotdir,plotfile)

stdout.write("Writing eps file: %s\n" % plotfile)
plt.write_eps(plotfile)
stdout.write("Converting to png\n")
dpi=90
esutil.misc.exec_process('converter -d %s %s' % (dpi,plotfile))

