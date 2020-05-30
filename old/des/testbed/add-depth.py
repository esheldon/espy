from __future__ import print_function
import numpy
import fitsio
import healpix_util as hu
import esutil as eu

def main():
    depth_file="/astro/u/esheldon/masks/des/sva1-gold/sva1_gold_1.0_nside4096-64_ring_i_weights.fits"

    depth_map=hu.readMap(depth_file)

    radec_file="sva1_coadd_spte-radec.fits"
    outfile=radec_file.replace('.fits','-depth.fits')

    print("reading:",radec_file)
    data=fitsio.read(radec_file)

    print("getting depth values")
    depth=depth_map.get_mapval(data['ra'], data['dec'])

    data=eu.numpy_util.add_fields(data, [('maglim','f4')])
    data['maglim'] = depth

    print("writing to:",outfile)
    fitsio.write(outfile, data, clobber=True)

main()
