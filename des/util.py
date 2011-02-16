"""
Processing SE run outputs:

    # create goodlist/badlist
    deswl.wlpipe.check_shear(serun, band)
    also use the /bin/check_shear.py

    # also generate the checkpsf pbs scripts with
    deswl.wlpipe.generate_se_pbsfiles(serun, band, type='checkpsf')
    



    #
    # QA plots
    #

    # look at averaged PSF patterns and star shape patterns
    des/bin/plot_checkpsf.py
    plot_checkpsf --serun wlse0003 --types allbothbin

    # Size mag diagrams
    des/bin/plot_size_mag.py


    # compare sizes from shapelet scale and the square dff
    # of sigma0 and sigma_psf
    des/bin/compare-sizes.py


    # All exposures from dc4
    des.util.dc4_plot_exposures(type='png')
    /bin/dc4_exposure_plot.py




    # Create sub-samples by region for dc4.
    des.util.dc4_create_shear_samples(serun, objclass, ftype=, delim=)
    /bin/dc4_create_samples.py


    

"""
import math
import os
import sys
from sys import stdout,stderr

import numpy
from numpy import array, where, arctan2
from esutil import wcsutil

import esutil
from esutil import json_util
from esutil import numpy_util
from esutil.plotting import setuplot, bwhiskers
from esutil.ostools import path_join, getenv_check


try:
    import deswl
    from deswl import wlpipe
except:
    stderr.write('Could not import deswl/wlpipe\n')


try:
    import pyfits
except:
    stderr.write('Could not import pyfits\n')


try:
    import recfile
except:
    pass


_defs={}
_defs['root'] = '/data/Archive'
_defs['project'] = 'DES'
_defs['ext'] = '.fits'
_defs['pixscale'] = 0.296 # arcsec/pixel



r2d = 180.0/math.pi
d2r = math.pi/180.0

element_sep = '-'
max_string_length = '125'
sx_number_width = 6
















bacon_ramin = 335.0
bacon_ramax = 340.0
bacon_decmin = -25.0
bacon_decmax = -20.0
    
def dc4_create_shear_samples(serun, objclass='gal',
                             ftype='rec',
                             delim=None, 
                             indir=None, 
                             outdir=None):
    """
    Read in the collated file and take sub-samples

    According to Huan this is the setup for DC4
    note ra=-25 is really ra=335.0

    These are in pixel coordinates
        x -> DEC
        y -> -RA
    But we measure in RA/DEC space.  

    RA          Dec        gamma1  gamma2
    -----------------------------------------
    >= -25 deg  >= -25 deg     0.05    0
    >= -25      <  -25        -0.025   0.025
    <  -25      >= -25         0.05    0.05
    <  -25      <  -25         0      -0.025

    rotated ra point is 335.0
    """

    rc=deswl.files.Runconfig(serun)
    header = rc.asdict()
    header['serun'] = header['run']
    del header['run']

    infile = deswl.files.wlse_collated_path(serun, 'all', ftype='rec', 
                                            dir=indir)


    stdout.write('Reading: %s\n' % infile)
    indata=esutil.io.read(infile)
    indata=indata.view(numpy.recarray)


    if ftype == 'fits':
        fitsh=pyfits.Header()
        fitsh.update('pyvers',header['pyvers'],'Version of python')
        fitsh.update('wlvers',header['wlvers'],'Version of wl pipeline')
        fitsh.update('tmvvers',header['tmvvers'],'Version of tmv')
        fitsh.update('esutilv',header['esutilvers'],'Version of esutil')
        fitsh.update('serun',header['serun'],'wl se run id')
        fitsh.update('dataset',header['dataset'],'data set used')




    # initial cuts
    stdout.write("Before shear cuts: %s\n" % indata.size)
    w0 = se_shear_cuts(indata, objclass)

    indata = indata[w0]

    stdout.write("After shear cuts: %s\n" % indata.size)




    stdout.write('Extracting regions\n')

    flist=[]
    regions=[1,2,3,4]
    for region in regions:
        stdout.write("\tregion: %s\n" % region)

        outfile=deswl.files.wlse_collated_path(serun,
                                               objclass,
                                               dir=outdir,
                                               region=region, 
                                               ftype=ftype, 
                                               delim=delim)

        w = dc4_region_select(indata.ra, indata.dec, region)
        stdout.write("\t\tFound: %s\n" % w.size)

        outdata = indata[w]

        stdout.write('\t\tWriting to file: %s\n' % outfile)

        if ftype == 'fits':
            # we're going to dump the data anyway, so we can noswap for
            # faster speed
            pyfits.writeto(outfile, outdata, header=fitsh, clobber=True)
        else:
            esutil.sfile.write(outdata, outfile, delim=delim, header=header)


    outfile=deswl.files.wlse_collated_path(serun,
                                           objclass,
                                           dir=outdir,
                                           region=regions, 
                                           ftype=ftype, 
                                           delim=delim)

    stdout.write('\tWriting combined file: %s\n' % outfile)
    if ftype == 'fits':
        pyfits.writeto(outfile, indata, header=fitsh, clobber=True)
    else:
        esutil.sfile.write(indata, outfile, delim=delim, header=header)


def collate_se_pointing(serun, exposurename):
    """
    Collate the output data from 
        $DESDATA/wlbnl/$SERUN/pointing
    With the input images and catalogs.  All I do is look at the .json
    files to get the images/catalogs and copy all into a single 
    directory under 
        $DESDATA/wlbnl/$SERUN/collated/pointings
    """

    import shutil

    dir = deswl.files.wlse_collated_dir(serun)
    dir = path_join(dir, 'pointings', exposurename)
    if not os.path.exists(dir):
        os.mkdir(dir)

    # now read all the json files and copy data over
    for ccd in range(1,62+1):
        statfile=deswl.files.wlse_path(exposurename, ccd, 'stat', serun=serun)

        stdout.write("\nReading stat file: %s\n" % statfile)
        stat = esutil.io.read(statfile)

        for type in ['imfile','catfile','psf','fitpsf','shear','stars']:
            fullpath=stat[type]

            if fullpath.find('scratch') != -1:
                fullpath = \
                    fullpath.replace(stat['rootdir'],deswl.files.des_rootdir())

            bname=os.path.basename(fullpath)
            outfile = path_join(dir, bname)

            stdout.write("Copying from %s to %s\n" % (fullpath,outfile))
            shutil.copy2(fullpath, outfile)
        

def plot_shearxy_byccd(serun, region, 
                       example_wcs_byccd=None, outfile=None, typ='pdf'):
    """
    Take the input data, split by ccd, and plot each separately.  The x,y
    are converted to a ra,dec type system that preserved the layout of
    the camera.  We don't use the actual ra/dec here because we want to
    be able to combine multiple exposures on the same plot.
    """

    from numpy import arctan2, pi

    fpath=deswl.files.wlse_collated_path(serun,'gal',ftype='rec',region=region)
    data = esutil.io.read(fpath, verbose=True)
    if data.size == 0:
        stdout.write("No objects in region: %s\n" % region)
        return


    outdir=os.path.dirname(fpath)
    outdir=os.path.join(outdir, 'plots')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile=os.path.basename(fpath).replace('.rec', '-shear-byccd.'+typ)
    outfile=os.path.join(outdir, outfile)
    stdout.write("Will write image: %s\n" % outfile)


    plt=setuplot('Agg')
    plt.clf()
    if example_wcs_byccd is None:
        example_wcs_byccd = get_exposure_wcs_example()

    h,rev = esutil.stat.histogram(data['ccd'], rev=True, min=1,max=62)

    allstats = []
    #for ccd in range(1,63):

    running_shear1=0.0
    running_shear2=0.0
    running_shear_weight=0.0
    for iccd in range(len(h)):
        if rev[iccd] != rev[iccd+1]:
            w = rev[ rev[iccd]:rev[iccd+1] ]

            ccd = data['ccd'][w[0]]

            stdout.write('\tccd=%s\n' % ccd)

            # no copy is made here!
            stats = stats_shearxy(data[w], nx=3, ny=3)

            ccd_wcs = example_wcs_byccd[ccd]
            xx,yy = ccd_wcs.image2sky(stats['mx'],stats['my'])

            xx -= 337.3
            yy += 15.0
            stats['mx'] = xx
            stats['my'] = yy

            allstats.append(stats)

            weights = 1.0/(stats['mshear1_err']**2 + stats['mshear2_err']**2)
            s1 = stats['mshear1']*weights
            s2 = stats['mshear2']*weights
            s1sum=s1.sum()
            s2sum=s2.sum()
            wsum = weights.sum()

            running_shear1 += s1sum
            running_shear2 += s2sum
            running_shear_weight += wsum

            thise1 = s1sum/wsum
            thise2 = s2sum/wsum
            thisangle = 0.5*arctan2(thise2, thise1)*180./pi
            stdout.write("\t\t<e1>=%s\n" % thise1 )
            stdout.write("\t\t<e2>=%s\n" % thise2 )
            stdout.write("\t\t<angle>=%s\n" % thisangle )


            #plt.plot(xx,yy,'.',markersize=1)

    #plt.show()

    stats = numpy_util.combine_arrlist(allstats)

    mshear1 = running_shear1/running_shear_weight
    mshear2 = running_shear2/running_shear_weight
    mangle = 0.5*arctan2(mshear2, mshear1)*180.0/pi
    mshear = numpy.sqrt( mshear1**2 + mshear2**2 )

    mangle_err = stats['mangle'].std()/numpy.sqrt(stats['mangle'].size)

    print "mins:",stats['mshear1'].min(), stats['mshear2'].min(), \
            stats['mshear'].min()
    print "maxs:",stats['mshear1'].max(), stats['mshear2'].max(), \
            stats['mshear'].max()
    stdout.write("overall averages: \n\tshear1: "
                 "%s\n\tshear2: %s\n\tangle: %s\n" % (mshear1, mshear2,mangle))
    stdout.write('\tangle error approximately %s\n' % mangle_err)

    # x component of the "vector" version
    u = stats['mshear']*numpy.cos(stats['mangle'])#*10000
    # y component of the "vector" version
    v = stats['mshear']*numpy.sin(stats['mangle'])#*10000

    # scale=1 means a vector of length 1 will cover essentially
    # all the plot.  I want them slightly smaller, so I'm using
    # scale=1.5
    scale=1.5
    ax=plt.axes()
    from matplotlib.ticker import MultipleLocator as ml
    ax.xaxis.set_minor_locator(ml(0.1))
    ax.yaxis.set_minor_locator(ml(0.1))
    #plt.quiver(stats['mx'], stats['my'], u, v,headwidth=0,
    #             scale=scale, pivot='middle')

    whiskers(plt, stats['mx'], stats['my'], u, v)


    xtext=336.0 - 337.3
    whiskers(plt, 336.18-337.3, -14.1+15.0, 0.05, 0.0, color='blue')
    plt.text(xtext, -14.1+15, "0.05", verticalalignment='center')


    ystart=-15.6 + 15
    ystep=-0.08
    istep = 0
    ytext=ystart+istep*ystep

    plt.text(xtext, ytext, r"$\langle \gamma_1 \rangle=%0.3f$" % mshear1,
               verticalalignment='center')
    istep+=1
    ytext=ystart+istep*ystep
    plt.text(xtext, ytext, r"$\langle \gamma_2 \rangle=%0.3f$" % mshear2,
               verticalalignment='center')
    istep+=1
    ytext=ystart+istep*ystep
    plt.text(xtext,ytext,r"$\langle \theta \rangle=%0.2f^{\circ}$" % mangle,
               verticalalignment='center')


    # plot a whisker representing the average
    istep+=1
    ytext=ystart+istep*ystep
    svec1 = mshear*numpy.cos( mangle*pi/180. )
    svec2 = mshear*numpy.sin( mangle*pi/180. )
    plt.text(xtext, ytext, r"$\langle \gamma \rangle=%0.3f$" % mshear, 
               verticalalignment='center')
    istep+=1
    ytext=ystart+istep*ystep
    whiskers(plt, 336.225-337.3, ytext, svec1, svec2, color='red')


    label = 'region%s' % region
    ax=plt.axes()
    plt.text(0.90, 0.9, label, 
               horizontalalignment='center', 
               verticalalignment='center', 
               transform=ax.transAxes, 
               fontsize=18)
    # this is so the arrows have the right angle

    plt.axis('equal')
    plt.ylim(-16.1+15.0,-13.9+15.0)

    if outfile is not None:
        stdout.write("Writing file: %s\n" % outfile)
        plt.savefig(outfile, bbox_inches='tight', pad_inches=0.2)
    else:
        plt.show()

    #print 'mean shear: ',stats['mshear']
    return example_wcs_byccd

def stats_xy_byccd(data, fields, nx=3, ny=3, typ='median', 
                   weights=None, index=None):

    if index is None:
        index=numpy.arange(data['x'].size, dtype='i4')

    ntot = nx*ny*62
    print 'total number of bins:',ntot

    dt=[('mx','f4'),('my','f4'),('ccd','i1')]
    for f in fields:
        if f is not 'x' and f is not 'y':
            new_dt = ('m'+f, 'f4')
            dt.append(new_dt)

    stats = numpy.zeros(ntot, dtype=dt)

    i=0
    for ccd in xrange(1,1+62):

        i1 = i*nx*ny
        i2 = i1+nx*ny

        stats['ccd'][i1:i2] = ccd
        w, = where(data['ccd'][index] == ccd)
        if w.size > 0:
            w=index[w]
            tstats = stats_xy(data, fields, nx=nx, ny=ny, typ=typ, 
                              weights=weights, index=w)

            for f in fields:
                #print f,i1,i2,stats['m'+f][i1:i2].shape,tstats['m'+f][:].shape
                stats['m'+f][i1:i2] = tstats['m'+f][:]
            if 'x' not in fields:
                stats['mx'][i1:i2] = tstats['mx'][:]
            if 'y' not in fields:
                stats['my'][i1:i2] = tstats['my'][:]

        i += 1
    return stats



def stats_xy(data, fields, nx=3, ny=3, typ='median', weights=None, index=None):

    if index is None:
        index=numpy.arange(len(data['x']), dtype='i4')

    dt=[('mx','f4'),('my','f4')]
    for f in fields:
        if f is not 'x' and f is not 'y':
            new_dt = ('m'+f, 'f4')
            dt.append(new_dt)

    h,rev = esutil.stat.histogram2d(data['x'][index], data['y'][index], 
                                    nx=nx, ny=ny, rev=True)

    stats = numpy.zeros(h.size, dtype=dt)

    for i in range(h.size):
        if rev[i] != rev[i+1]:

            # Note difference from IDL: no -1
            w = rev[ rev[i]:rev[i+1] ]

            # now w.r.t. original index
            w=index[w]

            for f in fields:
                fdata = data[f][w]

                if typ == 'mean':
                    if weights is not None:
                        m,e = esutil.stat.wmom(fdata,weights[w],calcerr=True) 
                    else:
                        m = fdata.mean()
                elif typ=='median': 
                    m=numpy.median(fdata)
                else:
                    raise ValueError("typ='mean' or 'median'")
                
                stats['m'+f][i] = m

            if 'x' not in fields:
                stats['mx'][i] = data['x'][w].mean()
            if 'y' not in fields:
                stats['my'][i] = data['y'][w].mean()

    return stats


def stats_shearxy(data, nx=20, ny=20):

    h,rev = esutil.stat.histogram2d(data['x'], data['y'], nx=nx, ny=ny, 
                                    rev=True)

    stats = numpy.zeros(h.size, dtype=[('mx','f4'),('my','f4'),
                                       ('mshear1','f4'),('mshear1_err','f4'),
                                       ('mshear2','f4'),('mshear2_err','f4'),
                                       ('mshear','f4'),('mangle','f4')])
    for i in range(h.size):
        if rev[i] != rev[i+1]:

            # Note difference from IDL: no -1
            w = rev[ rev[i]:rev[i+1] ]

            # errors seem unusually large, use a large SN softening
            SN=1.5
            e2 = ( SN**2 + data['shear_cov00'][w] + data['shear_cov11'][w] )
            weights = 1.0/e2
            m1, e1 = esutil.stat.wmom(data['shear1'][w],weights,calcerr=True)
            m2, e2 = esutil.stat.wmom(data['shear2'][w],weights,calcerr=True)

            meanx = data['x'][w].mean()
            meany = data['y'][w].mean()

            stats['mx'][i] = meanx
            stats['my'][i] = meany
            stats['mshear1'][i] = m1
            stats['mshear1_err'][i] = e1
            stats['mshear2'][i] = m2
            stats['mshear2_err'][i] = e2

            stats['mshear'][i] = numpy.sqrt( m1**2 + m2**2 )            
            stats['mangle'][i] = 0.5*numpy.arctan2(m2, m1)

    return stats


def dc4_plot_exposures(serun, overplot=True, prompt=False, typ='png', data=None):

    import astro_util
    plt = setuplot('Agg')

    plt.rcParams['figure.figsize'] = 8.5,11
    plt.clf()


    fpath = deswl.files.wlse_collated_path(serun, 'gal', ftype='rec', 
                                           region=[1,2,3,4])
    fields=['ra','dec','exposurename','ccd']
    if data is None:
        data=esutil.io.read(fpath, verbose=True, fields=fields, norecfile=True)

    if overplot:
        outdir=os.path.dirname(fpath)
        outdir=os.path.join(outdir, 'plots')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile=os.path.basename(fpath).replace('.rec', '-pos.'+typ)
        outfile=os.path.join(outdir, outfile)
        stdout.write("Will write image: %s\n" % outfile)


    stdout.write("shifting ra\n")
    sra=astro_util.shiftra(data['ra'])

    stdout.write("Getting ranges\n")
    minra=sra.min()
    maxra=sra.max()
    midra=(maxra+minra)/2.

    mindec=data['dec'].min()
    maxdec=data['dec'].max()
    middec=(maxdec+mindec)/2.


    from mpl_toolkits.basemap import Basemap
    #map = Basemap(projection='moll', lon_0=0)
    map = Basemap(projection='cass', 
                  llcrnrlon=-33, llcrnrlat=mindec,
                  urcrnrlon=-18, urcrnrlat=-5,
                  lon_0 = midra, lat_0=middec)


    #xrange=map(minra, maxra)
    #yrange=map(mindec,maxdec)


    #i=raw_input('hit a key: ')
    #plt.xlim(xrange[0], xrange[1])
    #plt.ylim(yrange[0], yrange[1])

    #stdout.write("Writing figure: %s\n" % outfile)
    #plt.savefig(outfile, dpi=100, bbox_inches='tight',pad_inches=0.5)
    #return

    stdout.write('Getting exposure names\n')
    enames = numpy_util.unique(data['exposurename'],values=True)

    stdout.write('Plotting exposures\n')

    ntot=len(enames)
    i=0
    for ename in sorted(enames):
        stdout.write('%s/%s %s\n' % (i,ntot,ename))
        i += 1
        w, = numpy.where( data['exposurename'] == ename )
        if w.size > 0:

            x, y = map(sra[w], data['dec'][w])
            if overplot:
                plt.plot( x, y, '.',markersize=2)
            else:
                plt.clf()
                plt.plot( x, y, '.',color='tan',markersize=2)

                ax=plt.axes()
                plt.text(0.85,0.9,ename,
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=ax.transAxes,fontsize=18)


                for ccd in range(1,63):
                    w2,=numpy.where(data['ccd'][w] == ccd)
                    if w2.size > 0:
                        w2=w[w2]

                        m=x.mean()
                        my=y.mean()

                        plt.text( mx, my, str(ccd), 
                                   fontsize=14, 
                                   weight='extra bold',
                                   horizontalalignment='center', 
                                   verticalalignment='center')



            if prompt:
                x=raw_input('hit a key: ')
                if x == 'q':
                    return



    #merid=numpy.arange(0.0, 365., 5.)
    #para =numpy.arange(-90., 95., 5.)
    merid=numpy.arange(0.0, 365., 5.)
    para =numpy.arange(-90., 95., 5.)

    merid_label=numpy.ones(merid.size)
    para_label=numpy.ones(para.size)
    stdout.write("drawing meridians\n")
    map.drawmeridians(merid, labels=[0,0,0,1], labelstyle='+/-')#, yoffset=-0.05)
    stdout.write("drawing parallels\n")
    map.drawparallels(para,labels=[1,0,0,0], labelstyle='+/-')#, xoffset=0.05)

    #plt.xlabel('RA')
    #plt.ylabel('DEC')

    if overplot:
        stdout.write("Writing figure: %s\n" % outfile)
        plt.savefig(outfile, 
                      dpi=100, 
                      transparent=False,
                      bbox_inches='tight',pad_inches=0.5)


def _get_exposure_wcs_byccd(run, exposurename, ccd):
    redfile=deswl.files.red_image_path(run, exposurename, ccd, check=True)
    stdout.write("Getting wcs from %s\n" % redfile)
    h = pyfits.getheader(redfile, ext=1)
    wcs = wcsutil.WCS(h)
    return wcs

def _get_redimage_wcs(redfile, verbose=False):
    if verbose:
        stdout.write("Getting wcs from %s\n" % redfile)
    h = pyfits.getheader(redfile, ext=1)
    wcs = wcsutil.WCS(h)
    return wcs


def get_exposure_wcs(exposurename, ccd=None, dataset=None, verbose=False):
    """
    wcs = get_exposure_wcs(run, exposurename, ccd=None)

    if ccd is None, a list for all ccds is returned.  If ccd is a scalar, 
    a single wcs object is returned
    """

    if dataset is None:
        dataset = deswl.wlpipe.default_dataset('wlse')

    band=deswl.wlpipe.getband_from_exposurename(exposurename)

    infodict = deswl.files.collated_redfiles_read(dataset, band)

    imdict={}
    for tinfo in infodict['flist']:
        if tinfo['exposurename'] == exposurename:
            tmp_ccd = int(tinfo['ccd'])
            imdict[tmp_ccd] = tinfo['imfile']


    if ccd is None:
        ccd = range(1,63)
    if not isinstance(ccd, list):
        return _get_redimage_wcs(imdict[ccd])

    all_wcs = {}
    for this_ccd in ccd:
        wcs = _get_redimage_wcs(imdict[this_ccd])
        wcs['ccd'] = this_ccd
        all_wcs[this_ccd] = wcs

    return all_wcs
        

def get_exposure_wcs_example():
    from esutil.ostools import getenv_check, path_join
    import cPickle as pickle

    exposurename='decam--24--15-i-6'
    d=getenv_check('DESFILES_DIR')
    pickle_path=path_join(d, 'wcsexample', 'wcs-'+exposurename+'.pickle')
    if not os.path.exists(pickle_path):
        stdout.write("Creating pickled version: %s\n" % pickle_path)
        tmpwcs = get_exposure_wcs(exposurename)
        pickle.dump(tmpwcs, open(pickle_path,'w'))
        stdout.write("Don't forget to svn add!\n")

    stdout.write("Loading wcs example: %s\n" % pickle_path)
    wcs = pickle.load(open(pickle_path))
    return wcs
































def ExtractPathInfo(name):
    """
    This makes assumptions about the directory structure
    """

    info = {}
    # remove the data dir
    data_dir = conf['data_dir']
    if data_dir[-1] != os.sep:
        data_dir += os.sep 
    reldir = name.replace(data_dir, '')


    if reldir == name:
        raise ValueError('File does not live in in tree %s' %
                         conf['data_dir'])

    relsplit = reldir.split(os.sep)
    desc = element_sep.join(relsplit)

    info['relpath'] = reldir
    info['filedesc'] = desc
    info['basename'] = os.path.basename(name)
    info['dirname'] = os.path.dirname(reldir)
    return info








def SXNumber2String(sx_number):
    return str(sx_number).rjust(sx_number_width, '0')

def ObjDesc(filedesc, number):
    numstr = SXNumber2String(number)
    objdesc = element_sep.join( [filedesc,numstr])
    return objdesc

def PGStuffFileInfo(fileclass, flist_filename, outdir=None):
    import npypg
    import es_util as eu

    # the output table is the files table
    files_table = conf['files_table']


    # info for this file class
    finfo = conf['files'][fileclass]

    # info about this table
    tabinfo = conf['tables'][files_table]

    # read the file list
    flist = open(flist_filename).readlines()
    nf = len(flist)

    # loop over the file names
    for fnum in range(nf):
        if True and (fnum > 0):
           return

        f = flist[fnum].strip()

        stdout.write( 'File (%d/%d): %s\n' % (fnum+1,nf,f))

        try:
            if not os.path.exists(f):
                raise RuntimeError('path does not exist: %s\n' % f)
        except:
            stderr.write('%s %s\n' % sys.exc_info()[0:2])
            continue

        # get some path information
        pathinfo = ExtractPathInfo(f)

        info = numpy.zeros(1, dtype=tabinfo['numpy_desc'])
        info['fileclass'] = finfo['fileclass']
        info['fileformat'] = finfo['fileformat']
        info['relpath'] = pathinfo['relpath']
        info['filedesc'] = pathinfo['filedesc']
        info['basename'] = pathinfo['basename']
        info['dirname'] = pathinfo['dirname']

        stdout.write('Stuffing file info\n') 
        try:
            npypg.array2table(info, files_table, 
                              unique=tabinfo['unique_cols'],
                              serial_extra=tabinfo['serial_cols'],
                              outdir=outdir)
        except:
            stderr.write('Error stuffing metadata for file %s\n' % f)
            stderr.write('%s %s\n' % sys.exc_info()[0:2])
            continue



def PGCatStuff(flist, stufftype, ftype='red_cat'):
    """
    Stuff the file info into the files table or the data into the obj table
    """


    import npypg
    import es_util as eu


    allf = open(flist).readlines()
    nf = len(allf)
    #for f in open(flist):
    fnum = 0

    finfo = conf["files"][ftype]
    tables = conf['tables']

    files_table = "files"

    for f in allf:
        if False and (fnum > 0):
           return
        f=f.strip()
        stdout.write( 'File (%d/%d): %s\n' % (fnum+1,nf,f))
        stdout.flush()

        try:
            if not os.path.exists(f):
                raise RuntimeError('path does not exist: %s\n' % f)
        except:
            stderr.write('%s %s\n' % sys.exc_info()[0:2])
            fnum+=1
            continue

        stdout.write( 'Reading data\n')
        try:
            if finfo['format'].lower() == 'fits':
                if finfo['type'] == 'bintable':
                    t,h = eu.fits2array(f, finfo['hdu'])
                elif finfo['type'] == 'header':
                    t = pyfits.getheader(f, finfo['hdu'])
                else:
                    raise ValueError("Don't support typoe: %s yet\n" % 
                                     finfo['type'])
            else:
                raise ValueError('Unsupported file format: %s\n' % finfo['format'])
        except:
            stderr.write('Error reading from file: %s\n' % f)
            stderr.write('%s %s\n' % sys.exc_info()[0:2])
            fnum+=1
            continue


        fnameinfo = ExtractPathInfo(f)

        if stufftype.lower() == 'files': 
            info = numpy.zeros(1, dtype=tables[files_table]['desc'])
            info['type'] = finfo['type']
            info['format'] = finfo['format']
            info['relpath'] = fnameinfo['relpath']
            info['desc'] = fnameinfo['desc']
            info['basename'] = fnameinfo['basename']
            info['dirname'] = fnameinfo['dirname']

            stdout.write('Stuffing metadata\n') 
            try:
                npypg.array2table(meta, files_table, 
                                  unique=['fileid','filedesc','filepath'], 
                                  serial_extra=['fileid'])
            except:
                stderr.write('Error stuffing metadata for file %s\n' % f)
                stderr.write('%s %s\n' % sys.exc_info()[0:2])
                fnum+=1
                continue
        elif stufftype.lower() == 'data':
            stdout.write( 'Stuffing cat data\n')
            try:
                q="SELECT fileid FROM %s WHERE filedesc = '%s'" % \
                        (files_table, finfo['desc'])
                stdout.write(q+'\n')
                res = npypg.query(q)
                sdef = 'S'+max_string_length
                #newdesc = [('filedesc',sdef),('objdesc',sdef)]
                newdesc = [('fileid','i8'),('objdesc',sdef)]
                new = eu.add_fields(t, newdesc) 

                new['fileid'] = res[0]['fileid']

                for i in range(len(new)):
                    new[i]['objdesc'] = ObjDesc(finfo['desc'], new[i]['NUMBER'])

                npypg.array2table(new, objtable, unique=['objid','objdesc'], 
                                  serial_extra=['objid'])
            except ValueError:
                stderr.write("Error stuffing file '%s' %s\n" % (f,err))
                stderr.write('%s %s\n' % sys.exc_info()[0:2])
                fnum += 1
                continue
        else:
            raise ValueError('Unsupported stuff type: %s\n' % stufftype)

        fnum+=1

    if stufftype.lower() == 'files':
        npypg.query('analyze %s' % files_table)
    elif stufftype.lower() == 'data':
        npypg.query('analyze %s' % objtable)
    stdout.write('\n')












def FitsHeaderAsDict(fname, ext=0):
    hdr=pyfits.getheader(fname, ext=ext)

    hdict={}
    for name,val in hdr.items():
        hdict[name.lower()] = val

    return hdict

        










































def ReadScampHead(fname, verbose=False):
    """
    Read one of these ascii .head files that SExtractor uses
    as a list of dictionaries
    """
    f=open(fname, "r")

    # Simpler to read the whole thing and then process it as needed
    data = f.readlines()
    f.close()

    # We will only process the keywords
    #   e.g.  x = y
    # and I won't worry about arrays
    
    headers=[]
    linedict={}
    for line in data:
        line=line.strip()
        # remove comments
        line = line.split('/')[0].strip()

        if line == 'END':
            headers.append(linedict)
            linedict={}
        else:
            # Now check for the = and use that to create an element in the
            # dictionary
            # for now only use lines with the = sign
            ls = line.split('=')
            if len(ls) == 2:
                # Execute the line but put a linedict['keyword'] around the
                # left hand side
                lhs = ls[0].strip().lower()
                lhs = "linedict['"+lhs+"']"
                rhs = ''.join( ls[1:] ).strip()

                comm = lhs+" = "+rhs
                if verbose:
                    sys.stdout.write('Executing: %s\n' % comm)
                try:
                    exec(comm)
                except:
                    if verbose:
                        sys.stdout.write('Command Failed: %s\n' % comm)
                        sys.stdout.write('  Stringifying\n')
                    rhs = "'" + rhs +"'"
                    comm = lhs+" = "+rhs
                    try:
                        if verbose:
                            sys.stdout.write('Executing: %s\n' % comm)
                        exec(comm)
                    except:
                        if verbose:
                            sys.stdout.write('Failed to read, skipping\n')
                            sys.stdout.write('   %s\n' % comm)


    return headers


def test_pv_invert_many_plot(rdict):
    plt = setuplot()
    plt.figure(1)
    plt.clf()
    plt.hist(rdict['rms'],20)
    plt.xlabel(r'Order increase')
    plt.ylabel('Number')
    plt.title('test')
    plt.show()

def test_pv_invert_many(find=True):
    import es_util
    from glob import glob
    sep = '-'*70

    dir=os.path.expanduser("~/data/astrometry")
    imdir=os.path.join(dir,'image')
    plotdir=os.path.join(dir,'plots')
    testdir=os.path.join(dir,'test')

    x = numpy.array([21.34, 1000.0, 1500.17], dtype='f8')
    y = numpy.array([11.21, 1000.0, 1113.92], dtype='f8')

    pattern = os.path.join(imdir,'*')
    images=glob(pattern)


    hdr=pyfits.getheader(images[0])
    n1=hdr['naxis1']
    n2=hdr['naxis2']
    xrang=numpy.array([1.0, n1], dtype='f8')
    yrang=numpy.array([1.0, n2], dtype='f8')
    #n=100
    n=10
    x,y = es_util.make_xy_grid(n, xrang, yrang)

    rms = []
    xdiffs = []
    ydiffs = []
    for imname in images:
        sys.stdout.write('image: %s\n' % imname)
        sys.stdout.write(sep+'\n')
        hdr=pyfits.getheader(imname)
        wcs=wcsutil.WCS(hdr)

        sys.stdout.write('doing xforms\n')
        ra,dec = wcs.image2sky(x,y)
        xn,yn = wcs.sky2image(ra,dec, find=find)

        xdiff = xn-x
        ydiff = yn-y

        t=(xdiff)**2 + (ydiff)**2
        trms = numpy.sqrt( t.sum()/t.size )
        rms.append(trms)
        xdiffs.append(xdiff)
        ydiffs.append(ydiff)

        sys.stdout.write('rms: %s\n' % trms)

    rms = numpy.array(rms,dtype='f8')
    out={}
    out['rms'] = rms
    out['x'] = x
    out['y'] = y
    out['xdifflist'] = xdiffs
    out['ydifflist'] = ydiffs
    return out


def test_invert_time():
    indir=os.path.expanduser("~/data/astrometry/image")
    plotdir=os.path.join(indir,'plots')
    imname=os.path.join(indir, "BCS2304-5436Og.070914_0304.088_04.fits")

    import timeit
    setup= ["import numpy",
            "import wcsutil",
            "import pyfits",
            "f='/Users/esheldon/data/astrometry/image/BCS2304-5436Og.070914_0304.088_04.fits'",
            "hdr=pyfits.getheader(f)",
            "wcs=wcsutil.WCS(hdr)",
            "x = numpy.arange(1000, dtype='f8')",
            "y = numpy.arange(1000, dtype='f8')",
            "ra,dec = wcs.image2sky(x,y)"]
    setup='; '.join(setup)

    n=3
    t=timeit.Timer("wcs.sky2image(ra,dec)",setup)
    sys.stdout.write("time for find: %s %s\n" % (t.timeit(n)/n))

    t=timeit.Timer("wcs.sky2image(ra,dec,find=False)",setup)
    sys.stdout.write("time for regular: %s %s\n" % (t.timeit(n)/n))


def test_invert(sip=False, invert=False, distort=True, 
                dofacrange=False, doorange=False):
    """
    Test the inversion routine in wcsutil
    """

    indir=os.path.expanduser("~/data/astrometry/image")
    plotdir=os.path.join(indir,'plots')
    if sip:
        imname=os.path.join(indir, 'newheader.fits')
        fend = '-sip'
    else:
        imname=os.path.join(indir, "BCS2304-5436Og.070914_0304.088_04.fits")
        fend = '--pv'
    
    catext=2

    sys.stdout.write("Image file: %s\n" % imname)

    hdr = pyfits.getheader(imname)

    wcs = wcsutil.WCS(hdr)

    if invert:
        if not dofacrange and not doorange:
            if sip:
                wcs.InvertSipDistortion()
            else:
                wcs.InvertPVDistortion()
        elif dofacrange:
            n = 20
            rms = numpy.zeros(n, dtype='f8')
            facs = list( range(1,n+1) )

            i=0
            for fac in facs:
                rms[i] = wcs.InvertDistortion(fac=fac,verbose=False)
            sys.stdout.write('%s\n' % facs)
            sys.stdout.write('%s\n' % rms)

            plt = setuplot()
            plt.clf()
            plt.plot(facs,rms,'.')
        elif doorange:
            fname='test-order-inverse'+fend+'.pdf'
            fname=os.path.join(plotdir,fname)
            n=6
            order_increases= list( range(n) )
            rms = numpy.zeros(n, dtype='f8')
            i=0
            for oi in order_increases:
                rms[i] = wcs.InvertDistortion(order_increase=oi,verbose=False)
                i += 1
            sys.stdout.write('order increases: %s\n' % order_increases)
            sys.stdout.write('rms: %s\n' % rms)

            plt = setuplot()
            plt.clf()
            plt.plot(order_increases, rms,'-.')
            #plt.xlabel(r'$\mathrm{order increase}$')
            #plt.ylabel('$\mathrm{rms}~[pix]$')
            #plt.rc('text',usetex=True)
            #plt.rc('font',**{'family':'sans-serif','sans-serif':['computer modern roman']})
            #plt.rc('font',**{'family':'serif','sans-serif':['Times New Roman']})

            from matplotlib.ticker import MultipleLocator as ml
            ax = plt.axes()
            ax.xaxis.set_minor_locator(ml(0.1))
            ax.yaxis.set_minor_locator(ml(0.1))

            plt.xlabel(r'Order increase')
            plt.ylabel(r'rms $[pix]$')
            ax.set_yscale('log')
            sys.stdout.write('Saving figure: %s\n' % fname)
            plt.savefig(fname)
            fname=fname.replace('.pdf','.eps')
            sys.stdout.write('Saving figure: %s\n' % fname)
            plt.savefig(fname)
            fname=fname.replace('.eps','.ps')
            sys.stdout.write('Saving figure: %s\n' % fname)
            plt.savefig(fname)

    else:
        x = numpy.array([21.34, 1000.0, 1500.17], dtype='f8')
        y = numpy.array([11.21, 1000.0, 1113.92], dtype='f8')

        #wcs.InvertDistortion(order_increase=1)
        ra,dec = wcs.image2sky(x,y)
        xp,yp = wcs.sky2image(ra,dec, find=True)
        sys.stdout.write('x,y= %s\n' % (x,y) )
        sys.stdout.write('crval= %s\n' % wcs.crval)
        sys.stdout.write('predicted ra,dec= %s\n' % (ra,dec))
        sys.stdout.write('inverted xp,yp= %s\n' % (xp,yp))
        sys.stdout.write('difference= %s\n' % (x-xp,y-yp))
        sys.stdout.write('fractional difference= %s\n' % ( (x-xp)/x,(y-yp)/y ) )




def test_image2sky_cat(local=True):
    """
    This one uses catalogs
    """

    if local:
        indir=os.path.expanduser("~/data/astrometry")
    else:
        indir="/data1/Archive/BCS/red/20080327000000_20070913/red/BCS2304-5436Og.070914_0304.088"
    imname=os.path.join(indir, "BCS2304-5436Og.070914_0304.088_04.fits")
    catname=os.path.join(indir,"BCS2304-5436Og.070914_0304.088_04_cat.fits")

    catext=2
    #wcs = FitsHeaderAsDict(imname, ext=catext-2)

    sys.stdout.write("Image file: %s\n" % imname)
    sys.stdout.write("Cat file: %s\n" % catname)

    wcs = pyfits.getheader(imname)
    cat = pyfits.getdata(catname, ext=catext).view(numpy.ndarray)


    wcsobj = wcsutil.WCS(wcs)

    diff = numpy.zeros(cat.size, dtype=[('ra','f8'),('radiff','f8'),
                                     ('dec','f8'),('decdiff','f8')])
    for objid in range(cat.size):
        ra = cat['ALPHA_J2000'][objid]
        dec = cat['DELTA_J2000'][objid]
        x = cat['X_IMAGE'][objid]
        y = cat['Y_IMAGE'][objid]

        era,edec = wcsobj.image2sky(x,y)
        sys.stdout.write("%s %s %s %s %s %s %s %s\n" % (x, y, ra, dec, era,
                                                        edec, (ra-era)*3600,
                                                        (dec-edec)*3600))

        diff['ra'][objid] = ra
        diff['radiff'][objid] = (ra-era)*3600.0
        diff['dec'][objid] = dec
        diff['decdiff'][objid] = (dec-edec)*3600.0

    sys.stdout.write('-'*70)
    sys.stdout.write('\nMean ra offset (arcsec): %s\n' %
                     (diff['radiff'].mean(),))
    sys.stdout.write('Mean abs(ra offset) (arcsec): %s\n' %
                     (numpy.abs(diff['radiff']).mean(),))
    sys.stdout.write('sdev ra offset (arcsec): %s\n' % (diff['radiff'].std(),))
    sys.stdout.write('Mean dec offset (arcsec): %s\n' %
                     (diff['decdiff'].mean(),))
    sys.stdout.write('Mean abs(dec offset) (arcsec): %s\n' %
                     (numpy.abs(diff['decdiff']).mean(),))
    sys.stdout.write('sdev dec offset (arcsec): %s\n' %
                     (diff['decdiff'].std(),))

    outfname=os.path.expanduser('~/tmp/test_wcs_cat.dat')
    sys.stdout.write('Writing to file: %s\n' % outfname)
    r=recfile.Open(outfname, "w", delim='\t')
    r.Write(diff)
    r.Close()





def test_john(method=1, imheader=False):
    ddir = os.path.expanduser('~/john_data')
    catfname = os.path.join(ddir, 'test.cat')
    imfname = os.path.join(ddir, 'test.fits')
    headfname = os.path.join(ddir, 'test.head')

    catext=2
    wcs = ReadScampHead(headfname)[catext-2]
    cat = pyfits.getdata(catfname, ext=catext).view(numpy.ndarray)
    imh = pyfits.getheader(imfname, catext-1)
    imwcs = {}
    for n,d in imh.items():
        imwcs[n.lower()] = d

    iobj = 20
    ra = cat['XWIN_WORLD'][iobj]
    dec = cat['YWIN_WORLD'][iobj]
    x = cat['XWIN_IMAGE'][iobj]
    y = cat['YWIN_IMAGE'][iobj]

    sys.stdout.write("Starting x,y = %s\n" % (x,y))

    if not imheader:
        wcss = wcsutil.WCS(wcs)
        nra,ndec = wcss.image2sky(x, y)
    else:
        wcss = wcsutil.WCS(imwcs)
        nra,ndec = wcss.image2sky(x, y)

    sys.stdout.write("%s %s %s %s %s %s\n" %
                     (ra,dec,nra,ndec,(ra-nra)*3600,(dec-ndec)*3600))









def GetWCS(image_id, verbose=False):
    query="""
        SELECT 
            wcsdim,
            ctype1,ctype2,
            cunit1,cunit2,
            crval1,crval2,
            crpix1,crpix2,
            cd1_1, cd1_2,
            cd2_1, cd2_2, 
            pv1_1, pv1_2, pv1_3, pv1_4, pv1_5, 
            pv1_6, pv1_7, pv1_8, pv1_9, pv1_10,
            pv2_1, pv2_2, pv2_3, pv2_4, pv2_5, 
            pv2_6, pv2_7, pv2_8, pv2_9, pv2_10
        FROM 
            image 
        WHERE 
            id=%d
        """ % image_id

    if verbose:
        sys.stdout.write('%s\n' % query)
    o=oracle_util.Connection()
    res = o.Execute(query)
    if res['pv1_3'] == -9999:
        res['pv1_3'] = 0.0
    if res['pv2_3'] == -9999:
        res['pv2_3'] = 0.0

    if 'pv1_0' not in res.dtype.names:
        res = AddFields(res, [('pv1_0','f8'),('pv2_0','f8')])
        res['pv1_0'] = 1.0
        res['pv2_0'] = 1.0

    return res


