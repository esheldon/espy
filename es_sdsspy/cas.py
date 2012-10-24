import os
import numpy

def make_fchart_pages(dir, ra, dec, radii, extra=None, nperpage=25):
    """
    parameters
    ----------
    dir:
        dir to hold files
    ra:
        array or scalar of ras
    dec:
        array or scalar of decs
    radii:
        array or scalar of radii
        in degrees
    """

    ra=numpy.array(ra, ndmin=1, dtype='f8')
    dec=numpy.array(dec, ndmin=1, dtype='f8')
    radii=numpy.array(radii, ndmin=1, dtype='f8')
    if radii.size == 1 and ra.size > 1:
        tmp = numpy.zeros(ra.size)
        tmp[:] = radii[0]
        radii=tmp

    dir=os.path.expandvars(dir)
    dir=os.path.expanduser(dir)
    if not os.path.exists(dir):
        os.makedirs(dir)

    npages=ra.size/nperpage
    left = ra.size % nperpage
    if left > 0:
        npages+=1

    indexf = os.path.join(dir,'index.html')
    print indexf
    index=open(indexf,'w')
    index.write("""
<html>
    <body>
""" )

    fold=None
    for i in xrange(npages):
        fname='index%03i.html' % i
        fname=os.path.join(dir,fname)
        base=os.path.basename(fname)
        index.write('<a href="{url}">{name}</a><br>\n'.format(url='./'+base,name=base))

        tr="""
            <tr><td>{ra}</td><td>{dec}</td><td><a href="{url}">url</a></td><td><a href="{urlnav}">nav</a></td></tr>"""

        topnav=''
        if i > 0:
            topnav='<a href="%s">previous</a> | ' % fold
        else:
            topnav += ''.join('&nbsp;'*8)
        if i < (npages-1):
            fnext='index%03i.html' % (i+1,)
            topnav += '<a href="%s">next</a>' % fnext

        
        print fname
        with open(fname,'w') as fobj:
            head="""
<html>
    <body>
        <p>
        {topnav}
        <p>
        <table border=1>
            <tr><th>RA</th><th>DEC</th><th>chart</th><th>navi</th></tr>
            """.format(topnav=topnav)
            fobj.write(head)
            i1 = i*nperpage
            i2 = (i+1)*nperpage
            tra=ra[i1:i2]
            tdec=dec[i1:i2]
            trad=radii[i1:i2]

            for ip in xrange(tra.size):
                scale=10.0 #arcsec/pixel
                radpix=trad[ip]*3600/scale # pixels
                width=int(2*radpix)
                url=fchart_url(tra[ip], tdec[ip], width=width, height=width,
                               scale=scale)
                urlnav=fchart_url(tra[ip], tdec[ip], width=width, height=width,
                                  scale=scale,nav=True)
                this_td=tr.format(ra=tra[ip],dec=tdec[ip],
                                  url=url,urlnav=urlnav)
                fobj.write(this_td)
            tail="""
        </table>
    </body>
</html>
            """
            fobj.write(tail)
            fobj.flush()
            fobj.close()
        fold=fname

    index.write("""
    </body>
</html>
                """)
    index.flush()
    index.close()
def fchart_url(ra, dec, width=400, height=400, opt='G', query='', 
               scale=0.4, nav=False):
    if nav:
        base_url='http://skyserver.sdss3.org/dr8/en/tools/chart/navi.asp'
    else:
        base_url="http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx"
        #?ra=154.266&dec=39.047&scale=0.79224&opt=&width=512&height=512
		#base_url = 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx'

    url = ("{base_url}?ra={ra}&dec={dec}&width={width}&height={height}"
           "&opt={opt}&scale={scale}&query={query}")
    url=url.format(base_url=base_url,
                   ra=ra,dec=dec,width=width,height=height,
                   opt=opt,scale=scale,
                   query=query)
    return url
