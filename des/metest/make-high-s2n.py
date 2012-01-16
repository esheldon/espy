"""
    %prog [options] run s2nmin

Description:

    get some high S/N objects and print out info to stdout, to be used for
    making cutouts.  Objects with shear_flags != 0 are removed.
"""
import os,sys
import des
import deswl

import esutil as eu
from esutil.numpy_util import where1
import numpy

from optparse import OptionParser

parser=OptionParser(__doc__)

options, args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

run = args[0]
s2nmin = float(args[1])

d=deswl.files.collated_dir(run)
cd=os.path.join(d,'cutouts')
if not os.path.exists(cd):
    os.makedirs(cd)

c=des.collate.open_columns(run)

flistfile=deswl.files.me_collated_path(run,'goodlist')
glist=eu.io.read(flistfile)

f=c['shear_flags'][:]
s2n = c['shear_s2n'][:]

logic = (f==0) & (s2n > s2nmin)

fi=c['input_flags'][:]
logic = logic & (fi==0)
fw=c['flags_weight'][:]
logic = logic & (fw==0)

w=where1(logic)

fields=['tilename','id','input_flags','shear_s2n','x_image','y_image']
data=c.read_columns(fields, rows=w)
fi=fi[w]

ind = numpy.array( ['%s-%05d' % (d['tilename'],d['id']) for d in data] )
s = ind.argsort()
data = data[s]


tile_old=None
fname=os.path.join(cd,'s2n-gt-%d.dat' % s2nmin)
print 'writing to:',fname
with open(fname,'w') as fobj:
    print >>fobj, 'tilename id input_flags x_image y_image shear_s2n image cat'
    for d in data:
        tile=d['tilename']
        if tile != tile_old:
            for g in glist:
                if g['tilename'] == tile:
                    break
        print >>fobj, d['tilename'],d['id'],d['input_flags'],d['x_image'],d['y_image'],d['shear_s2n'],g['image'],g['cat']
