"""
    %prog [options] run
"""
# note for SE we already cut on flags==0 I think
import os
import sys
import des
import deswl

import recfile

from esutil.numpy_util import where1

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-s","--stars", action='store_true', help="Use PSF stars")

options, args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run=args[0]

print 'run:',run
c=des.collate.open_columns(run)

d=deswl.files.collated_dir(run)

if options.stars:
    outfile='%s-stars-radec.dat' % run
    print 'getting star_flag'
    star_flag=c['star_flag'][:]
    print 'selecting stars'
    w=where1(star_flag != 0)
else:
    outfile='%s-radec.dat' % run

    print 'getting flags'
    flags   = c['shear_flags'][:]

    if run[0:2] == 'me':
        print 'doing flagsin and flags_weight for me'
        w=des.flags.select_good_me_bycol(c)
    else:
        w=des.flags.select_good_se_bycol(c)

print 'kept',w.size,'of',c['ra'].size
print 'reading ra/dec'

data=c.read_columns(['ra','dec'], rows=w)

outfile = os.path.join(d,outfile)
print 'writing:',outfile
with recfile.Recfile(outfile, mode='w', delim=' ') as fobj:
    fobj.write(data)
    
del data
del w


