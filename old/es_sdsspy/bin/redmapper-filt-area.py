"""
    %prog [options] mapfile infile outfile amin r1,r2,..,rn

Description
    Calculated the redmapper filtered area around the input ra,dec points.  
    
Inputs
    mapfile: 
        The STOMP map file
    infile:
        The input fits file.
            ra,dec: Equatorial coords
            da: Angular diameter distance to lens, for converting angles to
            megaparsecs.
    outfile:        
        The output fits file will be simply the filtered area fraction, matched row
        by row to the input file

    amin: 
        the minimum area to resolve in square degrees

        typical tycho star masks are about 6.e-4 square degrees.
        See the -p switch to alter how precisely this area is to be measured

    r1,r2,...,rn
        A comma-separated list of radii in the same units as da, e.g. Mpc
"""
import sys
from sys import stderr
import fitsio
import stomp
import numpy
from numpy import log10, sqrt, sin, cos
from numpy import pi as PI

import esutil as eu
from esutil.coords import randcap, sphdist

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--prec',default=0.1,
                  help="precision with which to determine amin. default: %default")
parser.add_option('-s','--show',action='store_true',
                  help="show a plot of random points")


class RandomCap:
    def __init__(self, amin, prec):
        """
        Generate random points on a spherical cap

        The density is such that the area of holes of size amin are measured
        with precision prec

        parameters
        ----------
        amin: minimum area to resolve in radians**2
        prec: The area of holes the size of amin are measured with this precision
        """
        self.amin=amin
        self.prec=prec
        self.density = 1/prec**2/amin

    def genrand(self, ra, dec, rad_deg):
        """
        Generate random points in the spherical cap

        parameters
        ----------

        ra,dec: scalars
            position in degrees
        rad:
            radius in degrees

        outputs:
            ra,dec,dist in degrees
        """

        rad = numpy.deg2rad(rad_deg)
        A = numpy.pi*rad**2
        nrand = self.get_nrand(A)

        ra_rand, dec_rand, dist = randcap(nrand,ra,dec,rad_deg, get_radius=True)
        return ra_rand, dec_rand, dist

    def get_nrand(self,A):
        """
        Get the number of randoms needed so that the area of holes of size
        self.amin are measured with precision self.prec
        """

        nrand = A*self.density
        return nrand

def show_and_wait(ra_cen,dec_cen,ra,dec,cont):
    import biggles
    plt=biggles.FramedPlot()

    pcen=biggles.Point(ra_cen,dec_cen,type='filled circle',size=2,color='magenta')
    plt.add(pcen)

    #p=biggles.Points(ra,dec,type='dot')
    p=biggles.Points(ra,dec,type='filled circle',size=0.5)
    plt.add(p)

    w,=numpy.where(cont == 0)
    if w.size > 0:
        pbad = biggles.Points(ra[w],dec[w],type='filled circle',size=0.5,color='red')
        plt.add(pbad)
        bt=biggles.PlotLabel(0.95,0.95,'Nmasked: %d' % w.size, halign='right')
        plt.add(bt)

    t=biggles.PlotLabel(0.1,0.95,'N: %d' % ra.size, halign='left')
    plt.add(t)

    plt.aspect_ratio=1
    plt.show()

    key=raw_input('hit a key (q to quit): ')
    if key.lower() == 'q':
        print >>stderr,'Exiting without writing'
        sys.exit(0)

def redmapper_filter(r_phys):
    # ESR put your radial filter here
    return numpy.ones(r_phys.size)

def process_points(map, data, rc, radii, show=False):
    import time
    dims = (len(radii),)
    output = numpy.zeros(data.size, dtype=[('fsum_tot','f8',dims),
                                           ('fsum_inmap','f8',dims),
                                           ('area_frac','f8',dims)])

    tm_rand=0.0
    tm_mask=0.0

    nobj=data.size
    for iobj,obj in enumerate(data):
        for irad,rad in enumerate(radii):

            if ((iobj+1) % 1000) == 0:
                print >>stderr,'%d/%d' % (iobj+1,nobj)

            rad_deg = numpy.rad2deg(rad/obj['da'])

            t0=time.time()
            ra_rand, dec_rand, dist = rc.genrand(obj['ra'],obj['dec'], rad_deg)
            r_phys = dist*obj['da']
            tm_rand += time.time()-t0

            t0=time.time()
            cont=map.Contains(ra_rand, dec_rand,"eq")
            tm_mask += time.time()-t0

            w,=numpy.where(cont==1)

            filt = redmapper_filter(r_phys)

            # ESR is this correct?
            fsum_tot   = filt.sum()
            fsum_inmap = filt[w].sum()

            if fsum_tot != 0:
                area_frac = fsum_inmap/fsum_tot
            else:
                area_frac = -9999

            output['area_frac'][iobj,irad] = area_frac
            output['fsum_tot'][iobj,irad] = fsum_tot
            output['fsum_inmap'][iobj,irad] = fsum_inmap

            if show:
                show_and_wait(obj['ra'],obj['dec'],ra_rand, dec_rand, cont)

    eu.misc.ptime(tm_rand,format='time for random generation: %s\n')
    eu.misc.ptime(tm_mask,format='time for mask checking: %s\n')
    return output

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 5:
        parser.print_help()
        sys.exit(45)

    mapfile=args[0]
    infile=args[1]
    outfile=args[2]
    amin=float(args[3])
    radii=args[4].split(',')
    radii=[float(r) for r in radii]
    prec=float(options.prec)
    show=options.show


    print >>stderr,'amin: %.3e' % amin
    print >>stderr,'precision:',prec
    print >>stderr,'radii:',radii
    print >>stderr,'outfile:',outfile

    print >>stderr,'reading mapfile:',mapfile
    map=stomp.Map(mapfile)

    print >>stderr,'reading infile:',infile
    data=fitsio.read(infile,rows=range(100))

    # ensure lower case due to mwrfits
    data.dtype.names = [n.lower() for n in data.dtype.names]

    amin_rad = amin*(PI/180.)**2
    rc = RandomCap(amin_rad,prec)
    print >>stderr,'density: %e per radian**2' % rc.density

    output=process_points(map, data, rc, radii, show=show)

    print >>stderr,'writing output:',outfile
    fitsio.write(outfile,output,clobber=True)

main()
