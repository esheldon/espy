import des
import sys
from sys import stdout

# "all" probably here
if len(sys.argv) < 2:
    stdout.write("usage: %s serun [objclass]\n" % sys.argv[0])
    stdout.write("  objclass defaults to 'all'\n")
    sys.exit(45)

serun=sys.argv[1]
objclass='all'
if len(sys.argv) > 2:
    objclass = sys.argv[1]

outdir='~/data/DES/wlbnl/%s/collated' % serun
stdout.write("Will write to local dir '%s' don't forget to "
             "copy to tutti\n" % outdir)
limit=None
des.util.collate_se_shear(serun, objclass, outdir=outdir, limit=limit)
