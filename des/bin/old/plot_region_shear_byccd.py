import des
import sys

serun='wlsetest0002'
typ=sys.argv[1]
wcs=None
for region in [1,2,3,4]:
    #for region in [1]:
    wcs=des.util.plot_shearxy_byccd(serun,region, example_wcs_byccd=wcs,
                                    typ=typ)

