import deswl
import sys

if len(sys.argv) < 3:
    sys.stderr.write("\nusage: check_shear.py serun band\n\n")
    sys.exit(1)

serun = sys.argv[1]
band = sys.argv[2]
deswl.wlpipe.check_shear(serun, band)
