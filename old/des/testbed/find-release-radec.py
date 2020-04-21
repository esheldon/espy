from __future__ import print_function
import numpy
import biggles
import fitsio
import desdb

def main():
    release='sva1_coadd_spte'
    conn=desdb.Connection()

    runs=desdb.files.get_release_runs(release)
    nrun=len(runs)

    withbands=["g","r","i","z"]
    nband=len(withbands)

    res=[]
    for i,run in enumerate(runs):
        if (i % 10) == 0:
            print("%d/%d" % (i+1,nrun))
        #print(run)

        query="""
        select
            tilename, band, ra, dec
        from
            coadd
        where
            run='%s'
        """ % run

        tres=conn.quick(query)

        rd=[r for r in tres if r['band'] in withbands]

        if len(rd) == nband:
            # they all have the same ra,dec
            res.append(rd[0])

    conn.close()

    ntile=len(res) 
    print("found",ntile,"tiles with bands",withbands)

    data=numpy.zeros(ntile, dtype=[('tilename','S13'),
                                   ('ra','f8'),
                                   ('dec','f8')])

    for i,r in enumerate(res):
        data['tilename'][i] = r['tilename']
        data['ra'][i]       = r['ra']
        data['dec'][i]      = r['dec']

    fitsname='%s-radec.fits' % release
    print(fitsname)
    fitsio.write(fitsname, data, clobber=True)

    plt=biggles.FramedPlot()
    pts=biggles.Points(data['ra'], data['dec'], type='filled circle', color='blue')

    plt.add(pts)
    plt.xlabel='RA'
    plt.ylabel='DEC'


    epsname='%s-radec.eps' % release
    print(epsname)
    plt.write_eps(epsname)

main()
