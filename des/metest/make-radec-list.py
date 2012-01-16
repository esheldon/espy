# note for SE we already cut on flags==0 I think
import os
import sys
import des
import deswl

import recfile

from esutil.numpy_util import where1

if len(sys.argv) < 2:
    print 'usage: make-radec-list run'
    sys.exit(1)

run=sys.argv[1]

print 'run:',run
d=deswl.files.collated_dir(run)
outfile = os.path.join(d,'%s-radec.dat' % run)

c=des.collate.open_columns(run)

print 'getting flags'
flags   = c['shear_flags'][:]

if run[0:2] == 'me':
    print 'doing flagsin and flags_weight for me'
    w=des.flags.select_good_me_bycol(c)
else:
    w=des.flags.select_good_se_bycol(c)

print 'kept',w.size,'of',flags.size
print 'reading ra/dec'
del flags

data=c.read_columns(['ra','dec'], rows=w)

print 'writing:',outfile
with recfile.Recfile(outfile, mode='w', delim=' ') as fobj:
    fobj.write(data)
    #pass
    
del data
del w


