import sys
import zphot

if len(sys.argv) < 2:
    print("usage: make-pofz-columns.py pzrun")
    sys.exit(45)

pzrun=sys.argv[1]

wc=zphot.weighting.WeightedOutputs()
wc.make_columns(pzrun)
