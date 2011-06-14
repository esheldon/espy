"""
    %prog pzrun

Description:
    Correct the p(z) using the ration N(z)/sum(p(z))

"""
from __future__ import print_function
import os
import sys
import zphot
import esutil as eu
from glob import glob

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    pzrun = args[0]

    corrfile = zphot.weighting.pofz_correction_file(pzrun)
    if not os.path.exists(corrfile):
        print("file not found:",corrfile)
        print("  generating...")
        zphot.weighting.make_pofz_correction(pzrun)


    print("reading correction:",corrfile)
    cstruct=eu.io.read(corrfile)

    pattern = zphot.weighting.pofz_file(pzrun, chunk='*')

    # this is just to get the chunk count
    files=glob(pattern)
    nchunk = len(files)
    for chunk in xrange(nchunk):
        data = zphot.weighting.read_pofz(pzrun, chunk)

        if cstruct.size != data['pofz'][0].size:
            raise ValueError("expected length %d for pofz but"
                             " got %d\n" % (cstruct.size, data['pofz'][0].size))

        data['pofz'] *= cstruct['corr']
        outfile=zphot.weighting.corrected_pofz_file(pzrun, chunk)
        d=os.path.dirname(outfile)
        if not os.path.exists(d):
            os.makedirs(d)

        print("Writing corrected p(z) to:",outfile)
        fobj=eu.recfile.Recfile(outfile, 'w', delim=' ')
        fobj.write(data)
        fobj.close()

main()
