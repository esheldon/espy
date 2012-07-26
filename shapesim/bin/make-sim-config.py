"""
    %prog [options] objmodel psfmodel

You must send shear
You must send psf_e1 and psf_e2 for psf models gauss and dgauss
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-e','--psf-ellip',default=None,
                  help="psf ellipticity e1,e2")
parser.add_option('-s','--shear',default=None,
                  help="shear g1,g2")
parser.add_option('--simname',default=None,
                  help="force sim name")
parser.add_option('--overwrite',action='store_true',
                  help="force overwrite")
parser.add_option('--nume',default='8',
                  help="number of e vals")
parser.add_option('--nums2',default='4',
                  help="number of s2 vals")

parser.add_option('--dryrun',action='store_true',
                  help="just print to the screen")

dg_template="""
name: %(simname)s

orient: ring
nring: 100

# we can trim later
dotrim: false

psfmodel: dgauss
psf_sigma: 1.4
psf_e1: %(psf_e1)s
psf_e2: %(psf_e2)s

psf_sigrat: 2.3
psf_cenrat: 0.09

objmodel: %(objmodel)s

shear: [%(shear1)s,%(shear2)s]

mins2: 0.5
maxs2: 2.0
nums2: %(nums2)s

mine:  0.05
maxe:  0.80
nume:  %(nume)s
"""

shortnames={'gauss':'g',
            'dgauss':'dg',
            'exp':'e',
            'dev':'d'}

def _get_simname(objmodel, psfmodel, i):
    simname = 'sim-%s%s%02d' % (shortnames[objmodel],shortnames[psfmodel],i)
    return simname
def get_simname_and_file(objmodel, psfmodel):

    i=1
    simname = _get_simname(objmodel, psfmodel, i)
    fname = shapesim.get_config_file(simname)
    while os.path.exists(fname):
        i+=1
        simname = _get_simname(objmodel, psfmodel, i)
        fname = shapesim.get_config_file(simname)

    return simname, fname

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    objmodel=args[0]
    psfmodel=args[1]

    if options.simname is not None:
        simname=options.simname
        fname = shapesim.get_config_file(simname)
        if os.path.exists(fname) and not options.overwrite:
            raise ValueError("file already exists: %s" % fname)
    else:
        simname, fname = get_simname_and_file(objmodel, psfmodel)

    if options.shear is None:
        parser.print_help()
        sys.exit(1)

    shear1,shear2=options.shear.split(',')
    if psfmodel in ['gauss','dgauss']:
        if options.psf_ellip is None:
            parser.print_help()
            sys.exit(1)
        psf_e1,psf_e2 = options.psf_ellip.split(',')
    else:
        raise ValueError("need to support more psf types")

    if objmodel == 'exp':
        text = dg_template % {'simname':simname,
                              'objmodel':objmodel,
                              'psf_e1':psf_e1,
                              'psf_e2':psf_e2,
                              'shear1':shear1,
                              'shear2':shear2,
                              'nums2':options.nums2,
                              'nume':options.nume}
    else:
        raise ValueError("need to support more obj types")

    print text
    print 'fname:',fname
    if not options.dryrun:
        with open(fname,'w') as fobj:
            fobj.write(text)
    else:
        print 'dryrun'
main()
