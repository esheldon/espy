"""
    %prog type sample

Description:

    Create an input catalog for objshear.  type must be 'scat' or 'lcat'.  The
    sample id implies a config file in ${ESPY_DIR}/lensing/config. 

    You can send -n to tell how many randoms to generate for random lcat
    catalogs.  Also required in that case is an extra name to add.

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-n",dest="nrand",default=None,
                  help="number of randoms to generate, default is from the config "
                       "you must send an extra name in this case")
parser.add_option("-e",dest="extra_name",default=None,
                  help="an extra name to add to the output file")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

type = args[0]
sample = args[1]
if type == 'scat':
    lensing.scat.create_input(sample)
elif type == 'lcat':
    if options.nrand is not None:
        lensing.lcat.create_input(sample, 
                                  nrand=int(options.nrand), 
                                  extra=options.extra_name)
    else:
        lensing.lcat.create_input(sample)
else:
    raise ValueError("type must be 'scat' or 'lcat'")

