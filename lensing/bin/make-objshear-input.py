"""
    %prog [options]

Description:

    Create an input catalog for objshear.  type must be 'scat' or 'lcat'.  The
    sample id implies a config file 

        ${ESPY_DIR}/lensing/config/lcat-{sample}.yaml. 

    You can send -n to tell how many randoms to generate for random lcat
    catalogs.  Also required in that case is an extra name to add.

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option("-t",dest="type",default=None,
                  help="type of input to make, required")
parser.add_option("-s",dest="sample",default=None,
                  help="sample name, required")

parser.add_option("-n",dest="nrand",default=None,
                  help="number of randoms to generate, default is from the config "
                       "you must send an extra name in this case")

parser.add_option("--split",default=None,
                  help=("split number 0 to nsplit-1."))

parser.add_option("-e",dest="extra_name",default=None,
                  help="an extra name to add to the output file")

options,args = parser.parse_args(sys.argv[1:])

if options.type is None or options.sample is None:
    parser.print_help()
    sys.exit(1)

type = options.type
sample = options.sample

split=options.split
if split is not None:
    split=int(split)

nrand=options.nrand
if nrand is not None:
    nrand=int(nrand)

if type == 'scat':
    # this will also run split() on the sample
    lensing.scat.create_input(sample=sample)
elif type == 'lcat':
    lensing.lcat.create_input(sample=sample, 
                              lens_split=split,
                              nrand=nrand,
                              extra=options.extra_name)
else:
    raise ValueError("type must be 'scat' or 'lcat'")

