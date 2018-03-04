from __future__ import print_function
import os
import numpy
import fitsio
import healpix_util as hu
import esutil as eu
import images

from scipy.interpolate import griddata

spte_ra_range=[60.,95.]
spte_dec_range=[-62.,-42.]

def get_map_radec():
    radec_file="map-with-radec-spte.fits"
    if not os.path.exists(radec_file):
        print("radec file",radec_file,"doesn't exist, creating it")

        #depth_file="/astro/u/esheldon/masks/des/sva1-gold/sva1_gold_1.0_nside4096-64_ring_i_weights.fits"
        dfile="/astro/u/esheldon/masks/des/neff/sva1_gold_1.0_nside4096-64_ring_i_neff.fits"
        print("reading:",dfile)
        dmap=hu.readDensityMap(dfile)
        print(dmap)

        print("npix:",dmap.data.size, dmap.hpix.npix)
        pixnums=numpy.arange(dmap.hpix.npix)

        print("getting ra,dec")
        ra,dec=dmap.hpix.pix2eq(pixnums)

        w,=numpy.where(  (ra > spte_ra_range[0])
                       & (ra < spte_ra_range[1])
                       & (dec > spte_dec_range[0])
                       & (dec < spte_dec_range[1]) )

        pixnums=pixnums[w]
        ra=ra[w]
        dec=dec[w]
        weight=dmap.data[w]
        weight *= (1./weight.max())

        dt=[('weight','f4'),('ra','f8'),('dec','f8')]
        data=numpy.zeros(ra.size, dtype=dt)
        data['ra']=ra
        data['dec']=dec
        data['weight']=weight

        print("writing:",radec_file)
        fitsio.write(radec_file, data, clobber=True)
    else:
        print("reading:",radec_file)
        data=fitsio.read(radec_file)

    return data
       
def get_map_image():
    data=get_map_radec()

    weight=data['weight']
    #minval=22.0
    #maxval=26.0
    #weight.clip(min=minval, max=maxval, out=weight)

    binsize=0.07
    res=eu.stat.histogram2d(data['ra'],
                            data['dec'],
                            z=weight,
                            #nx=200,
                            #ny=100)
                            xbin=binsize,
                            ybin=binsize)
    
    grid=res['zmean']
    grid_ra=res['xcenter']
    grid_dec=res['ycenter']

    return grid, grid_ra, grid_dec

def read_tiles():
    fname="sva1_coadd_spte-radec-depth.fits"
    print("reading:",fname)
    data=fitsio.read(fname)
    return data

testbed_tiles=[        'DES0511-5457','DES0516-5457',
               'DES0508-5540','DES0513-5540','DES0518-5540',
               'DES0509-5622','DES0514-5622','DES0520-5622',
                       'DES0511-5705','DES0517-5705',
              

               'DES0423-4748',
               'DES0453-4831',
               'DES0428-5914'
              ]
def main():
    import biggles
    biggles.configure('screen','width',1920)
    biggles.configure('screen','height',1200)

    tiledata=read_tiles()

    image, ra, dec = get_map_image()
    #images.view(image-image.min(),transpose=True)
    plt=images.view(image-image.min(),
                    xdr=[ra.min(), ra.max()],
                    ydr=[dec.min(), dec.max()],
                    transpose=False,
                    show=False)

    pts=biggles.Points(tiledata['ra'], tiledata['dec'],
                       type='filled circle', color='steel blue',
                       size=0.25)
    plt.add(pts)

    shoff=0.01
    #shoff=0.008
    for i in xrange(tiledata.size):
        ra=tiledata['ra'][i]
        dec=tiledata['dec'][i]
        tilename=tiledata['tilename'][i]
        if (    (spte_ra_range[0] < ra < spte_ra_range[1])
            and (spte_dec_range[0] < dec < spte_dec_range[1]) ):

            if tilename.strip() in testbed_tiles:
                color='green'
            else:
                color='red'
            dark_lab=biggles.DataLabel(ra,
                                       dec+0.05,
                                       tilename,
                                       color='black',
                                       size=0.3)

            plt.add(dark_lab)
            light_lab=biggles.DataLabel(ra-shoff,
                                        dec+0.05+shoff,
                                        tilename,
                                        color=color,
                                        size=0.3)
            plt.add(light_lab)

    plt.xrange=spte_ra_range
    plt.yrange=spte_dec_range
    plt.x1.ticks_style['color']='white'
    plt.x1.subticks_style['color']='white'
    plt.x1.spine_style['color']='white'

    plt.x2.ticks_style['color']='white'
    plt.x2.subticks_style['color']='white'
    plt.x2.spine_style['color']='white'

    plt.y1.ticks_style['color']='white'
    plt.y1.subticks_style['color']='white'
    plt.y1.spine_style['color']='white'

    plt.y2.ticks_style['color']='white'
    plt.y2.subticks_style['color']='white'
    plt.y2.spine_style['color']='white'

    epsname="map.eps"
    print("writing:",epsname)
    plt.write_eps(epsname)

    #pngname="map.png"
    #print("writing:",pngname)
    #plt.write_img(4096,2048,pngname)

main()
