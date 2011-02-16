def fchart_url(ra, dec, width=400, height=400, opt='G', query='', 
               nav=False):
    if nav:
        base_url = 'http://casjobs.sdss.org/dr7/en/tools/chart/navi.asp'
    else:
		base_url = 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx'

    url = ("{base_url}?ra={ra}&dec={dec}&width={width}&height={height}"
           "&opt={opt}&query={query}")
    url=url.format(base_url=base_url,
                   ra=ra,dec=dec,width=width,height=height,
                   opt=opt,query=query)
    return url
