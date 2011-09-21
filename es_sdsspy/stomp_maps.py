import sys
from sys import stdout
import os
import sdsspy
import numpy
from numpy import where, isscalar

try:
    import stomp
    INSIDE_MAP=stomp.cvar.Map_INSIDE_MAP
    FIRST_QUADRANT_OK=stomp.cvar.Map_FIRST_QUADRANT_OK
    SECOND_QUADRANT_OK=stomp.cvar.Map_SECOND_QUADRANT_OK
    THIRD_QUADRANT_OK=stomp.cvar.Map_THIRD_QUADRANT_OK
    FOURTH_QUADRANT_OK=stomp.cvar.Map_FOURTH_QUADRANT_OK

    QUAD1_OK = FIRST_QUADRANT_OK
    QUAD2_OK = SECOND_QUADRANT_OK
    QUAD3_OK = THIRD_QUADRANT_OK
    QUAD4_OK = FOURTH_QUADRANT_OK

    QUAD12_OK = QUAD1_OK + QUAD2_OK
    QUAD23_OK = QUAD2_OK + QUAD3_OK
    QUAD34_OK = QUAD3_OK + QUAD4_OK
    QUAD41_OK = QUAD4_OK + QUAD1_OK

    QUADALL_OK = QUAD1_OK + QUAD2_OK + QUAD3_OK + QUAD4_OK


except:
    pass

def quad_check(maskflags, strict=False, reverse=False):
    """
    if strict=False
        We will keep anything that has two adjacent quadrants contained
    if strict=True
        We demand *all* quandrants are unmasked
    """

    logic = quad_logic(maskflags, strict=strict, reverse=reverse)
    w,=where(logic)

    return w

def quad_logic(maskflags, strict=False, reverse=False):
    """
    if strict=False
        We will keep anything that has two adjacent quadrants contained
    if strict=True
        We demand *all* quandrants are unmasked
    """

    maskflags = numpy.array(maskflags, ndmin=1, copy=False)

    if not strict:
        good12 = (maskflags & QUAD12_OK) == QUAD12_OK
        good23 = (maskflags & QUAD23_OK) == QUAD23_OK
        good34 = (maskflags & QUAD34_OK) == QUAD34_OK
        good41 = (maskflags & QUAD41_OK) == QUAD41_OK

        if reverse:
            logic =  (good12==0) & (good23==0) & (good34==0) & (good41==0)
        else:
            logic =  good12 | good23 | good34 | good41

    else:
        if reverse:
            logic = (maskflags & QUADALL_OK) != QUADALL_OK
        else:
            logic = (maskflags & QUADALL_OK) == QUADALL_OK

    return logic

def quad_check_old(maskflags, strict=False):
    """
    if strict=False
        We will keep anything that has two adjacent quadrants contained
    if strict=True
        We demand *all* quandrants are unmasked
    """

    maskflags = numpy.array(maskflags, ndmin=1, copy=False)

    if not strict:
        check1 = \
            ( (maskflags & FIRST_QUADRANT_OK) != 0 ) & \
            ( (maskflags & SECOND_QUADRANT_OK) != 0 )
        check2 = \
            ( (maskflags & SECOND_QUADRANT_OK) != 0 ) & \
            ( (maskflags & THIRD_QUADRANT_OK) != 0 )
        check3 = \
            ( (maskflags & THIRD_QUADRANT_OK) != 0 ) & \
            ( (maskflags & FOURTH_QUADRANT_OK) != 0 )
        check4 = \
            ( (maskflags & FOURTH_QUADRANT_OK) != 0 ) & \
            ( (maskflags & FIRST_QUADRANT_OK) != 0 )
        w, = where( check1 | check2 | check3 | check4  )
    else:
        w, = where(  ( (maskflags & FIRST_QUADRANT_OK) != 0 )
                   & ( (maskflags & SECOND_QUADRANT_OK) != 0 )
                   & ( (maskflags & THIRD_QUADRANT_OK) != 0 )
                   & ( (maskflags & FOURTH_QUADRANT_OK) != 0 ) )

    return w




def map_dir():
    mapdir=os.getenv('MASK_DIR')
    if mapdir is None:
        raise ValueError("MASK_DIR is not defined")
    mapdir = os.path.join(mapdir, 'stomp-sdss')
    return mapdir

def map_name(name, maptype=None, band=None):
    mapdir=map_dir()

    fname=name

    if maptype is not None:
        fname += '-%s' % maptype
    if band is not None:
        fname += '-%s' % band

    fname += '.hmap'

    fname = os.path.join(mapdir, fname)
    return fname

def load(name, maptype=None, band=None, verbose=True):
    import stomp
    f=map_name(name, maptype=maptype, band=band)
    if not os.path.exists(f):
        raise IOError("stomp map file does not exist: %s" % f)
    stdout.write('Loading map: %s\n' % f)
    m = stomp.Map(f) 
    return m

def test(radius=None, n=1, map=None):
    import numpy
    import esutil
    if map is None:
        map=load('dr7','basic')

    ra = numpy.zeros(n, dtype='f8')
    dec = numpy.zeros(n, dtype='f8')
    ra[:] = 200
    dec[:] = 0

    maskflags = esutil.stomp_util.in_window(map, ra=ra, dec=dec, radius=radius)

    return map, maskflags


def read_boss_geometry():
    maskdir = map_dir()

    geom_file=os.path.join(maskdir,'boss_survey.par')

    if not os.path.exists(geom_file):
        raise ValueError("geometry file not found: '%s'" % geom_file)
    stdout.write("Reading geometry file: %s\n" % geom_file)
    bs = sdsspy.yanny.readone(geom_file)

    return bs

def plot_boss_geometry(color=None, colorwheel=None, plt=None, width=1, show=True,
                      region=None):
    """
    Plot the boundaries in the boss_survey.par file
    """
    import esutil as eu
    import biggles
    from biggles import FramedPlot, Curve

    bg = read_boss_geometry()

    if plt is None:
        plt = FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'

    if color is not None:
        colors = [color]*len(bg)
    elif colorwheel is not None:
        colors = colorwheel
    else:
        colors = ['red','blue','green','magenta','navyblue','seagreen',
                  'firebrick','cadetblue','green4']
        

    for i in xrange(len(bg)):
        b = bg[i]
        color = colors[i % len(colors)]
        c = eu.plotting.bbox( b['clambdaMin'], b['clambdaMax'], b['cetaMin'], b['cetaMax'],
                             color=color, width=width)
        plt.add(c)

    if region == 'ngc':
        plt.yrange = [-40.,50.]
        plt.xrange = [-80.,80.]
    elif region == 'sgc':
        plt.yrange = [105.,165.]
        plt.xrange = [-60.,60.]


    plt.aspect_ratio = (plt.yrange[1]-plt.yrange[0])/(plt.xrange[1]-plt.xrange[0])

    if show:
        plt.show()

    return plt

def create_boss_survey():
    """
    Create a stomp map from the boss_survey.par ceta,clambda bounds
    """
    import stomp

    maskdir = map_dir()

    map_file=map_name('boss','basic')
    stdout.write("Will write to file: %s\n" % map_file)

    geom_file=os.path.join(maskdir,'boss_survey.par')
    if not os.path.exists(geom_file):
        raise ValueError("geometry file not found: '%s'" % geom_file)
    stdout.write("Reading geometry file: %s\n" % geom_file)
    bs = sdsspy.yanny.read(geom_file,one=True)

    # convert all the bounds into stomp LatLon boundary objects,
    # then convert these to maps, and combine

    system = stomp.AngularCoordinate.Survey
    weight=1.0
    maxres=2048
    verbose=True
    mess = "-"*70 + "\nBound {i}/{ntot}: {lammin}:{lammax} {etamin}:{etamax}\n"
    for i in xrange(bs.size):
        lammin=bs['clambdaMin'][i]
        lammax=bs['clambdaMax'][i]
        etamin=bs['cetaMin'][i]
        etamax=bs['cetaMax'][i]
        stdout.write(mess.format(i=(i+1),
                                 ntot=bs.size,
                                 lammin=lammin,
                                 lammax=lammax,
                                 etamin=etamin,
                                 etamax=etamax))
        geom = stomp.LatLonBound(lammin,lammax,etamin,etamax,system)
        tmp_map = stomp.Map(geom,weight,maxres,verbose)

        if i == 0:
            map = tmp_map
        else:
            stdout.write("Ingesting new map\n")
            map.IngestMap(tmp_map)
        

    stdout.write("Writing map file: %s\n" % map_file)  
    map.Write(map_file)

def create_boss_survey_good():
    """
    Excise places where u amplifiers were not working
    """

    outfile = map_name("boss","good")
    stdout.write("Will write to file: %s\n" % outfile)

    map = load("boss", "basic")


	# ordered basically from bottom to top in eta

    ranges = numpy.zeros(5, dtype=[('lamrange','2f8'), ('etarange','2f8')])
    ranges['lamrange'][0] = [-66, -62.3]
    ranges['etarange'][0] = [-36.0, -35.6]

    ranges['lamrange'][1] = [-66, -30.5]
    ranges['etarange'][1] = [-35.6, -35.1]

    ranges['lamrange'][2] = [-18.8, -14.0]
    ranges['etarange'][2] = [-35.5, -35.3]

    ranges['lamrange'][3] = [-33.0, -30.5]
    ranges['etarange'][3] = [-30.8, -30.4]

    ranges['lamrange'][4] = [3.5,57.0]
    ranges['etarange'][4] = [-27.9, -27.55]

    system = stomp.AngularCoordinate.Survey
    weight=1.0
    maxres=2048
    verbose=True
    for i in xrange(ranges.size):

        stdout.write("lamrange: %s\n" % ranges['lamrange'][i])
        stdout.write("etarange: %s\n" % ranges['etarange'][i])
        bound = stomp.LatLonBound(ranges['lamrange'][i][0],
                                  ranges['lamrange'][i][1],
                                  ranges['etarange'][i][0],
                                  ranges['etarange'][i][1],
                                  system)


        badmap = stomp.Map(bound,weight,maxres,verbose)

        map.ExcludeMap(badmap)
        stdout.write('\n')

    stdout.write("Writing to file: %s\n" % outfile)
    map.Write(outfile)
	
def create_boss_survey_tycho():
    boss_tycho_file = map_name("boss","tycho")
    stdout.write("Will write to file: %s\n" % boss_tycho_file)

    stdout.write("Loading good\n")
    mgood = load("boss","good")
    stdout.write("Loading tycho\n")
    mtycho = load("tycho")

    stdout.write("Excluding tycho from good map\n")
    mgood.ExcludeMap(mtycho)

    stdout.write("Writing to file: %s\n" % boss_tycho_file)
    mgood.Write(boss_tycho_file)
