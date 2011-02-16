"""
Stuff relating to the tertiary archive
"""

import os
import esutil as eu
from esutil.ostools import path_join

def make_dc5b_data_html():
    d=os.environ['DESFILES_DIR']
    f = path_join(d,'dc5b','dc5b-runfiles-stat.json')

    runstat = eu.io.read(f)

    desdir = os.environ['DESDATA']

    www_dir=path_join(desdir,'www','dc5b')

    fhtml = path_join(www_dir, 'dc5bstat.html')
    fobj = open(fhtml, 'w')

    fobj.write('<html>\n')
    fobj.write('<body>\n')
    fobj.write('<head><link rel="stylesheet" href="../table.css" type="text/css"></head>\n')

    fobj.write('<h1>DC5b remote files</h1>\n')
    fobj.write("""<pre>
The following table summarizes the location of data for DC5b runs.
Only the desar machine is currently visible to BNL machines, so
only desar1,desar2,desardata locations included.
    </pre>""")

    fobj.write('<table class=simple>\n')
    fobj.write('    <tr><th>run</th><th>ftype</th><th>num</th><th>size (TB)</th><th>original<br>location</th><th>BNL<br>Has</th></tr>\n')

    # single epoch
    for ri in runstat['info']:
        run = ri['run']
        for ftype in ['red','red_cat']:
            num = ri[ftype]['num']
            sizetb = ri[ftype]['sizetb']
            host = ri[ftype]['host']
            have = ri[ftype]['bnlhas']

            fobj.write('    <tr><td>')

            data = [run, ftype, str(num), str(sizetb), host, have]
            data = '</td><td>'.join(data)
            fobj.write(data)

            fobj.write('</td></tr>\n')
    fobj.write('</table>\n')

    fobj.write('</html>\n')
    fobj.write('</body>\n')

    fobj.close()
