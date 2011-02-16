import zphot
import sys

files = sys.argv[1:]
zobj = zphot.cas.CasZphot()

for f in files:
    zobj.add_scinv_to_raw(f)
