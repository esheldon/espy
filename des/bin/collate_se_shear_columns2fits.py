import des
import sys
from sys import stdout

# "all" probably here
if len(sys.argv) < 2:
    stdout.write("usage: %s serun \n" % sys.argv[0])
    sys.exit(45)

serun=sys.argv[1]
des.collate.collate_se_shear_columns2fits(serun)
