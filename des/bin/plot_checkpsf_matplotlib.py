"""
    %prog [options] serun exposurename[s]

    example: %prog -d wlse0001 decam--25--10-i-1 decam--25--12-i-7 
"""
import sys
import os
from sys import stdout
from des import util as du
import deswl
import esutil
from esutil import numpy_util
from esutil.plotting import whiskers,polar2whisker,setuplot,set_minor_ticks

import numpy
from numpy import where

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-d",action="store_true",
        dest="diff",
        default=False,
        help="Generate a comparison plot and a diff plot")

parser.add_option("-n",dest="nbin",
                  default=3,
                  help="Number of x and y bins (nbinXnbin) Default %default")
parser.add_option("-f",dest="imageformat",
                  default='png',
                  help="Image file format. Default %default")

parser.add_option("-t",dest="ptypes",
    default='unbinned,comparebinned',
    help="Types of plots to make.  Can be a comma separated list. "
         "Can be currently be "
         "unbinned,compareunbinned,comparebinned,diffbinned.  Default %default")


def main():

    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    serun=args[0]
    # to put into a fiducial coordinate system
    exposurenames=args[1:]

    #
    diff = options.diff
    nbin=int(options.nbin)
    imageformat=options.imageformat
    ptypes = options.ptypes
    ptypes = ptypes.split(',')

    nx=nbin
    ny=nbin

    stdout.write("Reading example wcs\n")
    example_wcs_byccd = du.get_exposure_wcs_example()

    for exposurename in exposurenames:
        stdout.write('-'*30 + "\n")
        doplot(serun, exposurename, example_wcs_byccd, 
               imageformat=imageformat,
               nx=nx, ny=ny, ptypes=ptypes)
    


def doplot(serun, exposurename, example_wcs_byccd, nx=3,ny=3,diff=False,
           imageformat='png', ptypes=['unbinned','comparebinned']):

    stdout.write("serun: %s\n" % serun)
    stdout.write("exposurename: %s\n" % exposurename)
    stdout.write("\tReading checkpsf data\n")

    alldata = []
    allstats = []

    for ccd in range(1,62+1):

        try:
            data = deswl.files.wlse_read(exposurename,ccd,'checkpsf',
                                         serun=serun)

            fields=['e1','e2','e1interp','e2interp']
            stats = du.stats_xy(data, fields, nx=nx, ny=ny, typ='median')

            ccd_wcs = example_wcs_byccd[ccd]
            fx,fy = ccd_wcs.image2sky(data['x'], data['y'])
            mfx,mfy = ccd_wcs.image2sky(stats['mx'],stats['my'])

            newdata = numpy_util.add_fields(data, [('fx','f4'),('fy','f4')])
            newdata['fx'] = fx - 337.3
            newdata['fy'] = fy + 15.0

            fadd=[('mfx','f4'),('mfy','f4'),('me1diff','f4'),('me2diff','f4')]
            newstats = numpy_util.add_fields(stats,fadd)
            newstats['mfx'] = mfx - 337.3
            newstats['mfy'] = mfy + 15.0
            newstats['me1diff'] = newstats['me1interp']-newstats['me1']
            newstats['me2diff'] = newstats['me2interp']-newstats['me2']

            alldata.append(newdata)
            allstats.append(newstats)

        except:
            stdout.write("Failed ccd=%s\n" % ccd)
            print sys.exc_info()

    data = numpy_util.combine_arrlist(alldata)
    stats = numpy_util.combine_arrlist(allstats)

    stdout.write("\tPlotting whiskers\n")
    plt = du.setuplot('Agg')
    plt.clf()

    nplots = len(ptypes)
    xsize=8*nplots
    fig=plt.figure(figsize=(xsize,7))

    iplot = 1

    u, v = polar2whisker(data['e1'], data['e2'])
    uinterp, vinterp = polar2whisker(data['e1interp'], data['e2interp'])
    mu, mv = polar2whisker(stats['me1'],stats['me2'])
    muinterp, mvinterp = polar2whisker(stats['me1interp'],stats['me2interp'])
    mudiff, mvdiff = polar2whisker(stats['me1diff'],stats['me2diff'])



    # example size to plot
    psize=0.01

    if 'unbinned' in ptypes or 'compareunbinned' in ptypes:
        stdout.write("\t\tDoing plot type 'unbinned'\n")
        ax = fig.add_subplot(1,nplots,iplot)
        iplot += 1

        # plot individual stars
        scale=3.
        whiskers(ax, data['fx'], data['fy'], u, v, scale=scale,
                 linewidth=0.25)

        # a measure of scale
        xtext=-1.2
        ax.text(xtext, 0.9, str(psize), verticalalignment='center')
        whiskers(ax, xtext+0.18, 0.9, psize, 0.0, color='blue', scale=scale)

        if 'compareunbinned' in ptypes:
            stdout.write("\t\tOverplotting 'compareunbinned'\n")
            whiskers(ax, data['fx'], data['fy'], uinterp, vinterp, 
                     scale=scale, color='red', linewidth=0.25)
            
            # add a legend

            ypos = 1.1
            ystep = 0.07
            xtext = -1.1

            ax.text(xtext, ypos-ystep, "data", verticalalignment='center')
            whiskers(ax, xtext+0.23, ypos-ystep, psize, 0.0, 
                     color='black', scale=scale)
            ax.text(xtext, ypos-2*ystep, "interp", verticalalignment='center')
            whiskers(ax, xtext+0.23, ypos-2*ystep, psize, 0.0, 
                     color='red', scale=scale)


        ax.set_xlim(-1.3,1.3)
        ax.set_ylim(-1.3,1.3)
        set_minor_ticks(ax)

    if 'comparebinned' in ptypes:
        # comparison plot for the binned data
        ax = fig.add_subplot(1,nplots,iplot)
        iplot += 1

        scale=10.
        whiskers(ax, stats['mfx'], stats['mfy'], mu, mv, scale=scale)
        whiskers(ax, stats['mfx'], stats['mfy'], muinterp, mvinterp, 
                 scale=scale, color='red')

        # legend
        ypos = 1.1
        ystep = 0.07
        xtext = -1.1

        ax.text(xtext, ypos, str(psize), verticalalignment='center')
        whiskers(ax, xtext+0.23, ypos, psize, 0.0, color='blue', scale=scale)

        ax.text(xtext, ypos-ystep, "data", verticalalignment='center')
        whiskers(ax, xtext+0.23, ypos-ystep, psize, 0.0, 
                 color='black', scale=scale)
        ax.text(xtext, ypos-2*ystep, "interp", verticalalignment='center')
        whiskers(ax, xtext+0.23, ypos-2*ystep, psize, 0.0, 
                 color='red', scale=scale)

        ax.set_xlim(-1.3,1.3)
        ax.set_ylim(-1.3,1.3)
        set_minor_ticks(ax)

    if 'diffbinned' in ptypes:
        # comparison plot for the binned data
        ax = fig.add_subplot(1,nplots,iplot)
        iplot += 1

        scale=50.
        psize=0.001
        whiskers(ax, stats['mfx'], stats['mfy'], mudiff, mvdiff, scale=scale)

        ax.text(xtext, ypos, str(psize), verticalalignment='center')
        whiskers(ax, xtext+0.23, ypos, psize, 0.0, color='blue', scale=scale)

        ax.set_xlim(-1.3,1.3)
        ax.set_ylim(-1.3,1.3)
        set_minor_ticks(ax)


    outdir='/home/users/esheldon/www/tmp/plots'
    outdir=os.path.join(outdir, serun)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ptypes_string = '-'.join(ptypes)
    outfile = "%s-checkpsf-%s.%s" % (exposurename,ptypes_string, imageformat)
    outfile=os.path.join(outdir,outfile)
    stdout.write("\tWriting plot file: %s\n" % outfile)
    plt.savefig(outfile, bbox_inches='tight', pad_inches=0.1)

if __name__=='__main__':
    main()

