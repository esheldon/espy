"""
    %prog type sample

Description:
    Create an input catalog for objshear.  type must be 'scat' or 'lcat'.  The
    sample id implies a config file in ${ESPY_DIR}/lensing/config. 

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

type = args[0]
sample = args[1]
if type == 'scat':
    lensing.scat.create_input(sample)
elif type == 'lcat':
    lensing.lcat.create_input(sample)
else:
    raise ValueError("type must be 'scat' or 'lcat'")

