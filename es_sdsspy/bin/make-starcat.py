"""

Make a star catalog for clustering studies.

Trim the catalog to "intycho==1" and psfmag i (17.5,20.5)

I'm re-using flags for galaxies, which is not necessarily appropriate, but
since this doesn't need to be complete it shouldn't matter
"""

import os
import numpy
import es_sdsspy
import columns
import esutil as eu

#from optparse import OptionParser
#parser=OptionParser(__doc__)

MAGRANGE=(17.5,20.5)

def do_select(cols):
    select_fields=['psfflux_i','inbasic','objc_flags','objc_flags2','flags_r','flags_i']
    objs=cols.read_columns(select_fields, verbose=True)

    sel=es_sdsspy.select.Selector(objs)

    # already primary
    binned_logic = sel.binned_logic()
    object1_logic=sel.object1_logic()
    mask_logic=(objs['inbasic']==1)

    flux=objs['psfflux_i'].copy()
    flux.clip(0.001, flux.max(), flux)
    psfmag=22.5-2.5*numpy.log10(flux)
    mag_logic=(psfmag > MAGRANGE[0]) & (psfmag < MAGRANGE[1])


    logic=binned_logic & object1_logic & mag_logic & mask_logic

    w,=numpy.where(logic)
    print 'kept: %d/%d' % (w.size, objs.size)
    return w

def make_output(cols, w):
    d=os.path.expanduser('~/oh/starcat')
    f='dr8-starcat-basic-%.2f-%.2f.fits'
    f = f % MAGRANGE

    fname=os.path.join(d,f)

    dt=[('ra','f8'),
        ('dec','f8'),
        ('psfmag_i','f4'),
        ('ingood','i2'),
        ('instar','i2'),
        ('inbadfield','i2')]

    output=numpy.zeros(w.size, dtype=dt)

    flux=cols['psfflux_i'][w]
    flux.clip(0.001, flux.max(), flux)
    psfmag=22.5-2.5*numpy.log10(flux)

    output['psfmag_i'] = psfmag

    print 'copying'
    for n in output.dtype.names:
        if n != 'psfmag_i':
            print '    ',n
            output[n] = cols[n][w]

    print 'writing:',fname
    eu.io.write(fname,output,clobber=True)

def main():
    #options, args = parser.parse_args(sys.argv[1:])

    coldir=os.path.expanduser('~/sweep-reduce/dr8_final/primstar.cols')
    cols=columns.Columns(coldir)

    w=do_select(cols)
    make_output(cols, w)

main()
