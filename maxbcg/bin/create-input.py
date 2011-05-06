import sys
import maxbcg

if len(sys.argv) < 2:
    sys.stderr.write('usage: create-input.py procrun [filetype]\n')
    sys.stderr.write("    e.g. procrun='prim03' filetype='fits'\n")
    sys.stderr.write("    if filetype is not sent, both 'fits' and "
                     "'rec' are written\n")
    sys.exit(45)

procrun = sys.argv[1]

mcs = maxbcg.select.Selector(procrun)

mcs.select()
if len(sys.argv) > 2:
    filetype = sys.argv[2]
    mcs.write(filetype)
else:
    mcs.write('fits')
    mcs.write('rec')
