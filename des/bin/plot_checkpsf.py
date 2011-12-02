"""
    %prog [options] exposurename[s]

"""
import sys
import os
from sys import stdout,stderr
from des import util as du
import deswl
import esutil
from esutil import numpy_util
from esutil.ostools import path_join
from esutil.plotting import bwhiskers,polar2whisker

import numpy
from numpy import where
from numpy import sqrt

import biggles

PSF_VALIDATION = 0x10000000


from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-d",action="store_true",
        dest="diff",
        default=False,
        help="Generate a comparison plot and a diff plot")

parser.add_option("--serun",
                  default=None,
                  help="Use the given serun")

parser.add_option("--coldir",
                  default=None,
                  help="Use the input columns dir instead of the default")


parser.add_option("--dir",
                  default=None,
                  help="Run on individual files in the given dir")


parser.add_option("--magrange",
                  default=None,
                  help="limit the psf stars to this mag range. "
                       "example: --magrange=17,18. Default %default")


parser.add_option("--screen",action="store_true",
        dest="screen",
        default=False,
        help="Send plots to the screen additioin to a file. Default %default")

parser.add_option("-n",dest="nbin",
                  default=3,
                  help="Number of x and y bins (nbinXnbin) Default %default")

parser.add_option("--dpi",dest="dpi",
                  default=100,
                  help="DPI of png files. Default %default")

parser.add_option("-f",dest="imageformat",
                  default='png',
                  help="Image file format. Default %default")

parser.add_option("--scale",dest="scale",
                  default=5.0,
                  help="Scale the whiskers by this number %default")

parser.add_option("--types",dest="ptypes",
    default='starnobin,psfnobin,bothnobin,starbin,psfbin,bothbin',
    help="Types of plots to make.  Can be a comma separated list. "
         "If allstarbin, allbothbin, or allpsfbin, all exposures are used.  "
        "Can be currently be "
         "starnobin,psfnobin,bothnobin,starbin,psfbin,allbothbin,allstarbin,allpsfbin  Default %default")



def run_from_dir():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    exposurenames=args[:]
    dir=options.dir

    nbin=int(options.nbin)
    ptypes = options.ptypes
    ptypes = ptypes.split(',')

    screen = options.screen
    scale = options.scale
    dpi = options.dpi


    plotdir_base = path_join(dir, 'plots/psfcheck')
    if not os.path.exists(plotdir_base):
        os.makedirs(plotdir_base)


    dirs = {}
    allowed_types = \
        ['starnobin','psfnobin','bothnobin','starbin','psfbin','bothbin']


    for type in allowed_types:
        plotdir=path_join(plotdir_base, type)
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        dirs[type] = plotdir


    stdout.write("Reading example wcs\n")
    ewcs = du.get_exposure_wcs_example()

    for exposurename in exposurenames:

        stdout.write('-'*30 + "\n")

        stdout.write('%s\n' % exposurename)

        data = read_exposure_data(dir,exposurename)


        # this is an index into the psf stars
        wpsf, = where(data['psf_flags'] == 0)
        if wpsf.size == 0:
            raise ValueError("No objects passed psf_flags")
        else:

            fields = ['e1','e2','e1interp','e2interp']
            stats = du.stats_xy_byccd(data,fields, nx=nbin,ny=nbin, index=wpsf)

            for type in allowed_types:
                if type in ptypes:
                    # read these from the main table

                    if 'nobin' in type:
                        sx=data['x'][wpsf]
                        sy=data['y'][wpsf]
                        se1=data['e1'][wpsf]
                        se2=data['e2'][wpsf]
                        se1i=data['e1interp'][wpsf]
                        se2i=data['e2interp'][wpsf]
                        sccd=data['ccd'][wpsf]
                        sscale=scale
                    else:
                        sx=stats['mx']
                        sy=stats['my']
                        se1=stats['me1']
                        se2=stats['me2']
                        se1i=stats['me1interp']
                        se2i=stats['me2interp']
                        sccd=stats['ccd']
                        sscale=3*scale


                    plotdir=dirs[type]
                    psfile="%s-checkpsf-%s.eps" % (exposurename,type)
                    psfile=path_join(plotdir, psfile)
                    if 'star' in type:
                        plt=doplot(sx, sy, 
                                   se1,se2,
                                   sccd, ewcs, scale=sscale)
                    elif 'psf' in type:
                        plt=doplot(sx, sy, 
                                   se1i,se2i,
                                   sccd, ewcs, scale=sscale)
                    elif 'both' in type:
                        plt=doplot(sx, sy, 
                                   se1,se2,
                                   sccd, ewcs, 
                                   scale=sscale)
                        plt=doplot(sx, sy, 
                                   se1i,se2i,
                                   sccd, ewcs, 
                                   scale=sscale, 
                                   color='red',
                                   plt=plt)
                        # add the legend.  Only easy way to do this is
                        # to add fake ones off screen
                        xx=numpy.array([-9999.9,-9999.0],dtype='f8')
                        yy=numpy.array([-9999.9,-9999.0],dtype='f8')
                        cstars = biggles.Curve(xx,yy, color='black',linewidth=2)
                        cstars.label = 'stars'
                        cpsf = biggles.Curve(xx,yy, color='red',linewidth=2)
                        cpsf.label = 'psf'
                        key = biggles.PlotKey(0.1,0.9,[cstars,cpsf])
                        plt.add(key)

                    if screen:
                        print 'showing in a window'
                        plt.show()

                    stdout.write("Writing eps file: %s\n" % psfile)
                    plt.write_eps(psfile)
                    stdout.write("Converting to png\n")
                    esutil.misc.exec_process('converter -d %s %s' % (dpi,psfile))



def read_psf_stars_by_exposurename(cols, exposurename):
    # remember uid is the same as the row number for the main table
    if exposurename=='all':
        # we will use *all* psf stars
        stdout.write("Reading from all exposures\n")
        psf_flags = cols['psfstars']['psf_flags'][:]

        # we choose validation stars
        wpsf,=numpy.where( (psf_flags & PSF_VALIDATION) != 0)
        psf_uid = cols['psfstars']['uid'][wpsf]
        return psf_uid, wpsf
    else:
        all_uid = (cols['exposurename'] == exposurename)

        if all_uid.size == 0:
            raise ValueError("No exposure name '%s' found\n" % exposurename)
        else:
            # extract the psf objects that are a subset of these ids
            psf_uid, wpsf = cols['psfstars']['uid'].match(all_uid, both=True)

            # get the psf flags for these objects
            psf_flags = cols['psfstars']['psf_flags'][wpsf]

            # now get the ones with psf_flags set to zero
            psf_keep, = numpy.where(psf_flags == 0)
            if psf_keep.size == 0:
                stdout.write("No objects passed psf_flags cut")
                wpsf=numpy.array([],dtype='i4')
                psf_uid=numpy.array([],dtype='i4')
            else:
                wpsf = wpsf[psf_keep]
                psf_uid = psf_uid[psf_keep]

            return psf_uid, wpsf

    
def doplot(x, y, e1, e2, ccd, ewcs, scale=5.0, 
           plt=None, color='black'):
    """
    Just plot the e1/e2 whiskers for the input objects
    """

    if plt is None:
        plt = biggles.FramedPlot()

    ra = x.copy()
    dec = x.copy()
    ra[:] = -9999
    dec[:] = -9999
    for this_ccd in xrange(1,62+1):
        # put the x,y from each ccd onto the camera location
        # of our fiducial pointing, held in ewcs

        ccd_wcs = ewcs[this_ccd]

        wccd, = where(ccd == this_ccd)
        if wccd.size > 0:
            tra, tdec = ccd_wcs.image2sky(x[wccd],y[wccd])

            tra -= 337.3
            tdec += 15.0
            ra[wccd] = numpy.float32( tra )
            dec[wccd] = numpy.float32( tdec)

    # ok, let's make a plot

    # add whiskers for each e1/e2 to the plot
    u,v = polar2whisker(e1,e2)
    plt = bwhiskers(ra, dec, u, v, plt=plt, color=color, scale=scale)

    one_percent = 0.01*scale

    one_percent_x = numpy.array([0.9,0.9],dtype='f4')
    one_percent_y = numpy.array([1.0,1.0],dtype='f4')
    one_percent_x[0] -= one_percent*0.5
    one_percent_x[1] += one_percent*0.5
    one_percent_curve = biggles.Curve(one_percent_x,one_percent_y,color='blue',
                                      linewidth=2)
    plt.add(one_percent_curve)

    one_percent_label = biggles.PlotLabel(0.935,0.885,"0.01")
    plt.add(one_percent_label)

    plt.xrange = -1.3,1.3
    plt.yrange = -1.3,1.3
    plt.xlabel = r'$\Delta$RA'
    plt.ylabel = r'$\Delta$DEC'

    return plt



def run_from_serun():

    import columns
    options, args = parser.parse_args(sys.argv[1:])

    serun=options.serun

    nbin=int(options.nbin)
    ptypes = options.ptypes
    ptypes = ptypes.split(',')

    screen = options.screen
    scale = float(options.scale)
    dpi = options.dpi

    magrange = options.magrange
    if magrange is not None:
        magrange = numpy.fromstring(magrange,sep=',')
        if magrange.size != 2:
            raise ValueError("magrange should have two elements")


    allowed_types = \
        ['starnobin','psfnobin','bothnobin','starbin','psfbin','bothbin',
         'allstarbin','allpsfbin','allbothbin']

    # special case where we use all exposures
    if 'allstarbin' in ptypes or 'allpsfbin' in ptypes or 'allbothbin' in ptypes:
        exposurenames = ['all']
    else:
        if len(args) < 1:
            parser.print_help()
            sys.exit(45)

        exposurenames=args[:]




    # we will use this example set of wcs to put things onto
    # a camera plane

    stdout.write("Reading example wcs\n")
    ewcs = du.get_exposure_wcs_example()

    collated_dir = deswl.files.wlse_collated_dir(serun)
    plotdir_base = path_join(collated_dir, 'plots/psfcheck')
    if not os.path.exists(plotdir_base):
        os.makedirs(plotdir_base)
    dirs = {}

    for type in allowed_types:
        plotdir=path_join(plotdir_base, type)
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        dirs[type] = plotdir


    if options.coldir is not None:
        coldir=options.coldir
    else:
        coldir = deswl.files.coldir(serun)

    # open the column database
    stdout.write("Opening coldir: '%s'\n" % coldir)
    cols = columns.Columns(coldir)

    for exposurename in exposurenames:

        stdout.write('-'*30 + "\n")

        stdout.write('%s\n' % exposurename)

        # get index into main table for objects with this exposurename
        # note, this is equivalent to the uid column

        psf_uid, wpsf = read_psf_stars_by_exposurename(cols,exposurename)

        
        if psf_uid.size > 0:

            if magrange is not None:
                stdout.write("Limiting to mag range: %s" % magrange)
                imag = cols['imag'][psf_uid]
                w,=numpy.where( (imag > magrange[0]) & (imag < magrange[1]) )
                if w.size == 0:
                    raise ValueError("No objects in mag range: %s" % magrange)
                psf_uid = psf_uid[w]
                wpsf=wpsf[w]

            stdout.write("Found %s objects\n" % psf_uid.size)
            # ok, we now have all the indices we need to do our work
            # since psf_uid is an index back into the main table
            # and wpsf is row-by-row the corresponding index into the
            # psfstars table

            stdout.write("Reading x,y,ccd,psfe1,psfe2 from main table\n")
            x = cols['x'][psf_uid]
            y = cols['y'][psf_uid]
            ccd = cols['ccd'][psf_uid]
            
            e1interp = cols['interp_psf_e1'][psf_uid]
            e2interp = cols['interp_psf_e2'][psf_uid]

            # reading psfstars info
            stdout.write("Reading psf star e1,e2\n")
            psfstars_e1 = cols['psfstars']['e1'][wpsf]
            psfstars_e2 = cols['psfstars']['e2'][wpsf]

            stdout.write("Getting Stats\n")
            tdata={}
            tdata['x'] = x
            tdata['y'] = y
            tdata['ccd'] = ccd
            tdata['e1'] = psfstars_e1
            tdata['e2'] = psfstars_e2
            tdata['e1interp'] = e1interp
            tdata['e2interp'] = e2interp
            fields = [f for f in tdata.keys() if f != 'ccd']
            stats = du.stats_xy_byccd(tdata,fields,
                                      nx=nbin,ny=nbin)

            for type in allowed_types:
                if type in ptypes:
                    # read these from the main table

                    if 'nobin' in type:
                        sx=x
                        sy=y
                        se1=psfstars_e1
                        se2=psfstars_e2
                        se1i=e1interp
                        se2i=e2interp
                        sccd=ccd
                    else:
                        sx=stats['mx']
                        sy=stats['my']
                        se1=stats['me1']
                        se2=stats['me2']
                        se1i=stats['me1interp']
                        se2i=stats['me2interp']
                        sccd=stats['ccd']


                    plotdir=dirs[type]
                    psfile="%s-checkpsf-%s.eps" % (exposurename,type)
                    if magrange is not None:
                        psfile=psfile.replace('.eps','-mag%0.1f-%0.1f.eps' % magrange)
                    psfile=path_join(plotdir, psfile)
                    if 'star' in type:
                        plt=doplot(sx, sy, 
                                   se1,se2,
                                   sccd, ewcs, scale=scale)
                    elif 'psf' in type:
                        plt=doplot(sx, sy, 
                                   se1i,se2i,
                                   sccd, ewcs, scale=scale)
                    elif 'both' in type:
                        plt=doplot(sx, sy, 
                                   se1,se2,
                                   sccd, ewcs, 
                                   scale=scale)
                        plt=doplot(sx, sy, 
                                   se1i,se2i,
                                   sccd, ewcs, 
                                   scale=scale, 
                                   color='red',
                                   plt=plt)
                        xx=numpy.array([-9999.9,-9999.0],dtype='f8')
                        yy=numpy.array([-9999.9,-9999.0],dtype='f8')
                        cstars = biggles.Curve(xx,yy, color='black',linewidth=2)
                        cstars.label = 'stars'
                        cpsf = biggles.Curve(xx,yy, color='red',linewidth=2)
                        cpsf.label = 'psf'
                        key = biggles.PlotKey(0.1,0.9,[cstars,cpsf])
                        plt.add(key)

                    if screen:
                        print 'showing in a window'
                        plt.show()

                    stdout.write("Writing eps file: %s\n" % psfile)
                    plt.write_eps(psfile)
                    stdout.write("Converting to png\n")
                    esutil.misc.exec_process('converter -d %s %s' % (dpi,psfile))





def read_exposure_data(dir,exposurename):

    out_dtype=[('ccd','i1'),
               ('x','f4'),('y','f4'),
               ('psf_flags','i4'),
               ('e1','f4'),('e2','f4'),
               ('e1interp','f4'),('e2interp','f4')]

    datalist=[]
    for ccd in xrange(1,1+62):
        psf_file = path_join(dir,'%s_%02d_psf.fits' % (exposurename,ccd))
        psf=esutil.io.read(psf_file)
        psfe1 = psf['shapelets'][:,3]*sqrt(2)
        psfe2 = -psf['shapelets'][:,4]*sqrt(2)

        shear_file = path_join(dir,'%s_%02d_shear.fits' % (exposurename,ccd))
        shear=esutil.io.read(shear_file)
        e1interp = shear['interp_psf_coeffs'][:,3]*sqrt(2)
        e2interp = -shear['interp_psf_coeffs'][:,4]*sqrt(2)

        mpsf,mshear=numpy_util.match(psf['id'], shear['id'])

        tdata = numpy.zeros(mpsf.size, dtype=out_dtype)

        tdata['ccd'] = ccd
        tdata['x'] = psf['x'][mpsf]
        tdata['y'] = psf['y'][mpsf]
        tdata['e1'] = psfe1[mpsf]
        tdata['e2'] = psfe2[mpsf]

        tdata['psf_flags'] = psf['psf_flags'][mpsf]

        tdata['e1interp'] = e1interp[mshear]
        tdata['e2interp'] = e2interp[mshear]

        datalist.append(tdata)


    data = numpy_util.combine_arrlist(datalist)
    return data



def main():
    options, args = parser.parse_args(sys.argv[1:])

    if options.dir is not None:
        run_from_dir()
    elif options.serun is not None:
        run_from_serun()
    else:
        stderr.write("You must send --dir or --serun\n")
        parser.print_help()
        sys.exit(45)



if __name__=='__main__':
    main()

