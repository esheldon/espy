"""
Plot average yi values in bins of x

x could be shear_s2n, y1 y2 could be shear1 shear2
"""
import os,sys
import numpy

from esutil.numpy_util import where1
import esutil as eu
import columns
import converter
import des

args=sys.argv[1:]
if len(args) < 2:
    print 'usage: shearplot run x y1 y2 ...'
    sys.exit(1)


run = sys.argv[1]
fields = sys.argv[2:]

c=des.collate.open_columns(run)
print 'columns dir:',c.dir

pdir='$DESDATA/wlbnl/%s/collated/shear-vs' % run
pdir=os.path.expandvars(pdir)

use_input_flags=True
use_flags_weight=True

print 'reading shear_flags'
f=c['shear_flags'][:]

name_end=[]
logic = (f == 0)
if use_input_flags:
    print 'reading input_flags'
    fi=c['input_flags'][:]
    logic = logic & (fi == 0)
    name_end += ['iflags']

if use_flags_weight:
    print 'reading flags_weight'
    fw=c['flags_weight'][:]
    logic = logic & (fw == 0)
    name_end += ['wflags']

name_end = '-'.join(name_end)
print 'checking logic'
w=where1(logic)
print 'keeping:',w.size,'of',f.size

print 'reading columns:',fields
data=c.read_columns(fields, rows=w)

#s2n=c['shear_s2n'][w]
#err=numpy.sqrt(2)/s2n
#weights=1.0/(err**2 + 0.32**2)

weights=None
if fields[0] == 'shear_s2n':
    xlog=True
else:
    xlog=False

clip=False

keys={}
if 'se' in run and fields[0] == 'imag':
    keys['max'] = 14.5

print 'doing bhist_vs'
plts = eu.plotting.bhist_vs(data, *fields, nperbin=100000, xlog=xlog, 
                            weights=weights, clip=clip, **keys)

print 'writing out plots'
xfield=fields[0]
for field,plt in zip(fields[1:],plts):
    epsname='%s-%s-%s-%s' % (run,xfield,field,name_end)
    epsname=os.path.join(pdir,epsname)
    eu.ostools.makedirs_fromfile(epsname)
    if weights is not None:
        epsname+='_w'
    elif clip:
        epsname += '_clip'

    epsname+='.eps'
    print epsname
    plt.write_eps(epsname)
    converter.convert(epsname, verbose=True)
