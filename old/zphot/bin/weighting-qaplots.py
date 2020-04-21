"""
    %prog wrun

If neither -z or -v are sent then /both/ are performed.
"""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-s","--subind",default=None,
                  help="Do a subset of the indices for quick testing [0,subind)")
parser.add_option("-z","--zhist",default=None,action='store_true',
                  help="make the z histogram")
parser.add_option("--pzrun",default=None,
                  help="The pzrun for overplotting summed individual p(z)")
parser.add_option("-v","--varhist",default=None,action='store_true',
                  help="make the var comparision histogram")


options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

wrun = args[0]

subind = options.subind
if subind is not None:
    subind = int(subind)
    subind = range(0,subind)

pzrun = options.pzrun

do_zhist = options.zhist
do_varhist = options.varhist

if do_zhist is None and do_varhist is None:
    do_zhist=True
    do_varhist=True

cw = zphot.weighting.CompareWeighting(wrun, pzrun=pzrun)

if do_zhist:
    cw.zhist(dopng=True)

if do_varhist:
    cw.varhist(subphoto=subind)
