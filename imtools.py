#!/usr/bin/env python

import os
import sys
import re

from PIL import Image

image_pattern = 'jpg|jpeg|png|gif'
image_re = re.compile('.*('+image_pattern+')$',re.IGNORECASE)

_space='&nbsp;'

def GetImageFileList(dir):
    tdir = os.path.expanduser(dir)
    tdir = os.path.expandvars(tdir)
    files = os.listdir(tdir)

    image_files=[]
    for f in files:
        if image_re.match(f):
            f = os.path.join(dir,f)
            image_files.append(f)
    image_files.sort()
    return image_files

def ThumbDir(dir):
    return os.path.join(dir,'.thumbnails')

def ResizeName(fullpath,ftype='thumb'):
    dir, file = os.path.split(fullpath)
    front, ext = os.path.splitext(file)
    newdir = ThumbDir(dir)

    if ftype == 'thumb':
        f = front+'_thumb'+ext
    elif ftype == 'reduced':
        f = front+'_reduced'+ext
    else:
        sys.stdout.write("Unknown ftype: %s\n" % ftype)
        sys.exit(45)

    f = os.path.join(newdir, f)
    return f

def ResizeNames(flist, ftype='thumb', dodirs=False):
    newflist=[]
    for f in flist:
        newf = ResizeName(f,ftype=ftype)
        newflist.append(newf)
    return newflist

def ExtractSize(size_string):
    try:
        sx,sy = size_string.lower().split('x')
        sx,sy = int(sx),int(sy)
    except:
        sys.stdout.write("Size string is incorrectly formatted: %s\n" %
                         size_string)
    return (sx,sy)

def MakeThumbnailDir(dir):

    td = ThumbDir(dir)
    try:
        if not os.path.exists(td):
            os.mkdir(td)
    except:
        sys.stdout.write('Failed to create directory: %s\n' % td)
        sys.exit(45)


def ResizeDirImages(dir, size, ftype='thumb', overwrite=False):
    flist = GetImageFileList(dir)
    MakeThumbnailDir(dir)
    MakeResizedImages(flist, size, ftype=ftype, overwrite=overwrite)

def MakeResizedImages(flist, size, ftype='thumb', overwrite=False):
    """
    ResizeImage(flist, size, ftype)
    ftype: 'thumb' or 'reduced'

    A resized copy is placed in the .images2html directory if the original
    size is greater than the new size, else a symlink is put in
    """

    newflist = ResizeNames(flist,ftype=ftype)
    ntotpix = size[0]+size[1]
    for f,newf in zip(flist,newflist):
        mess='   Creating: '

        skip=False
        if os.path.exists(newf):
            if not overwrite:
                sys.stdout.write('   Reduced image already exists: %s\n' % newf)
                skip=True
            else:
                mess='   Overwriting: '
                os.remove(newf)
        if not skip:
            image = Image.open(f)
            osx, osy = image.size
            totpix = osx+osy
            if totpix > ntotpix:
                image.thumbnail( size, Image.ANTIALIAS ) 
                sys.stdout.write('%s %s\n' % (mess,newf))
                image.save(newf)
            else:
                mess = '   Image is smaller than requested size. Linking'
                sys.stdout.write('%s %s\n' % (mess,newf))
                os.symlink(f, newf)

def index_name(dir, page, npage, znum):
    if page == 0:
        name='index.html'
    elif page < 0:
        return ""
    elif page > (npage-1):
        return ""
    else:
        name='index'+str(page).zfill(znum)+'.html'

    name=os.path.join(dir, name)
    return name

def MakeWebpages(conf):
    """
    Loops over files in the entered config, which must contain
    flist, thumblist, redlist at the minimum
    """
    dir=conf['dir']
    imlist=conf['imlist']
    tlist=conf['thumblist']
    rlist=conf['redlist']

    tdir = './.thumbnails'

    nrow=conf['nrow']
    ncol=conf['ncol']
    nperpage = nrow*ncol
    nimage = len(imlist)
    npage = nimage/nperpage
    nleft = nimage % nperpage

    if nleft > 0 and npage != 1:
        npage+=1

    sys.stdout.write("# per page: %s\n" % nperpage)
    sys.stdout.write("# images: %s\n" % nimage)
    sys.stdout.write("# pages: %s\n" % npage)
    sys.stdout.write("# remainder: %s\n" % nleft)

    indent='  '
    znum=2
    if npage > 99:
        znum=3
    iim=0


    firstIndex = index_name(dir, 0, npage, znum)
    lastIndex  = index_name(dir, npage-1, npage, znum)
    firstIndex = '<a href="%s">First</a>' % firstIndex
    lastIndex = '<a href="%s">Last</a>' % lastIndex

    for page in xrange(npage):

        prevIndex = index_name(dir, page-1, npage, znum)
        iname=index_name(dir, page, npage, znum)
        nextIndex = index_name(dir, page+1, npage, znum)

        ifobj=open(iname,'w')

        ifobj.write("<html>\n")
        ifobj.write("<!-- Created by images2html Erin Sheldon -->\n")
        ifobj.write("<head>\n")
        ifobj.write("<link rel=\"STYLESHEET\" type=\"text/css\" href='"+conf['cssfile']+"'>\n")
        ifobj.write("</head>\n")

        ifobj.write("<body>\n")

        ifobj.write('<div id="content">\n')
        ifobj.write(indent+'<div class="thumbtable">\n')
        for row in range(nrow): 
            ifobj.write(indent*2+'<div class="row">\n')
            for col in range(ncol):
                ifobj.write(indent*3+'<div class="td">\n')
                crap,image = os.path.split(imlist[iim])
                crap,thumb = os.path.split(tlist[iim])
                crap,red = os.path.split(rlist[iim])

                turl = os.path.join(tdir,thumb)
                rurl = os.path.join(tdir,red)

                ifobj.write(indent*4+'<img src="'+turl+'">\n')

                ifobj.write(indent*4+'<br><span class="thumbname">'+image+'</span>\n')
                #ifobj.write(indent*4+'<font color=white>'+image+'</font>\n')

                ifobj.write(indent*3+'</div>\n')

                iim+=1

                # break out of column loop
                if iim > (nimage-1):
                    break

            # row end div
            ifobj.write(indent*2+'</div>\n')

            # break out of row loop
            if iim > (nimage-1):
                break

        # The "previous" and "next" links
        if nextIndex != "":
            nextIndex='<a href="%s">Next</a>' % nextIndex
        else:
            nextIndex=_space*4
        if prevIndex != "":
            prevIndex='<a href="%s">Previous</a>' % prevIndex
        else:
            prevIndex=_space*8



        ifobj.write(indent+'</div>\n')

        if page == 0:
            firstlast = lastIndex
        elif page == (npage-1):
            firstlast = firstIndex
        else:
            firstlast = firstIndex + '/' + lastIndex

        prevnextdiv="""
  <div class="prevNext">
    <div class="prevNextRow">
      <div class="prevtd">
        {prevIndex}
      </div>
      <div class="firstlast">
        {firstlast}
      </div>
      <div class="nexttd">
        {nextIndex}
      </div>
    </div>
  </div>\n""".format(prevIndex=prevIndex, firstlast=firstlast, nextIndex=nextIndex)

        ifobj.write(prevnextdiv)


        ifobj.write('</div>\n')
        ifobj.write('</body>\n</html>\n')
        ifobj.close() 





stylesheet="""
    /* Colors used: 
        #595B30  A kind of brown for the main text color
        #F0E68C  Khaki for the main background
        #898B60  A tannish color for nav background a few others
        #ffe4c4  bisque. visited links in navigation bar
    */

    /******************************************/
    /* Styling rules for the body             */
    /******************************************/

    html,body 
    {
      margin:0;

      color:#595B30;      /* A kind of brown */
      background:#000000; /* khaki; */
    }

    a:link {color:#595B30; text-decoration:none}
    a:visited {color:#595B30; text-decoration:none}
    a:hover {text-decoration:underline}

    /*-----------------------------------------------*/
    /* The main content, with absolute positioning   */
    /*-----------------------------------------------*/

    /*
    #content 
    {
      position:absolute;
      top:3em;
      left:8em;
    }
    */

    /*-----------------------------------------------*/
    /* A navigation bar will be on the left          */
    /*-----------------------------------------------*/

    #navigation 
    {
      width:8em;
      height:100%;
      padding:0.5em;

      background:#000000; /* A tannish color */
    }

    #navigation hr { display:none }
    #navigation h2 { color:#595B30 }
    #navigation a:link {color:white}
    #navigation a:visited {color:#ffe4c4} /* bisque */


    /*--------------------------------------------------------*/
    /* A class for the main table division on the index pages */
    /*--------------------------------------------------------*/

    .thumbtable 
    {
      display:table; 
      border-collapse:separate;

      margin:10px auto;


      width:640px;
      height:480px;

      /* this puts a background, and the border spacing, combined with 
         a background color of the td's, creates the border. I switched to
         putting an actual border on, see the row div.td */
      /*border-spacing:2px;*/
      /*background:#898B60;*/
    }

    .thumbname 
    { 
      font-style:italic;
      font-size:small;

      color:#898B60;
    }

    /*-----------------------------------------------*/
    /* rows and data classes                         */
    /*-----------------------------------------------*/

    .row 
    {
      display:table-row;
    }

    .row div 
    {
      display:table-cell;
      border-spacing:2px;
      height:33%; 

      background:#000000;

    }

    /* td is a subclass I guess */

    .row div.td 
    {
      text-align:center;
      vertical-align:middle;
      width:25%;

      /* these replace the effect of border-spacing and background in
         the main thumbtable. Better, since \"missing\" images are filled
         with the body background */
      border: 1px #898B60 solid;
    }

    .row div.td img { border-color:#595B30 }

    /*-----------------------------------------------*/
    /* A class for the bottom \"previous\" and \"next\"  */
    /* index                                         */
    /*-----------------------------------------------*/

    .prevNext 
    {  
      margin-left:auto;
      margin-right:auto;
      margin-top:0px;
      margin-bottom:0px;
      width:640px;
      /*border: 2px black solid;*/
      /*height:5em;*/
    }

    .prevNextRow 
    {  
      display:table-row;
    }

    .prevNextRow div
    { 
      display:table-cell;
    }

    .prevNextRow div.prevtd
    {  
      /*border: red 2px solid;*/
      width:1%;
      text-align:left;
      font-style:italic;
    }
    .prevNextRow div.nexttd
    {  
      /*border: red 2px solid;*/
      width:1%;
      text-align:right;
      font-style:italic;
    }

    .prevNextRow div.firstlast
    {  
      /*border: red 2px solid;*/
      width:100%;
      text-align:center;
      font-style:italic;
    }





    /*-----------------------------------------------*/
    /* A class for containing the images             */
    /*-----------------------------------------------*/

    .imagetable 
    {
      display:table; 
      border-collapse:separate;

      margin:10px auto;
      border-spacing:2px;
    }

    .imagetable img { border-color:#898B60 }

"""

def WriteCSS(file):
    f=open(file,'w')
    f.write(stylesheet)
    f.close()



_imhtml="""<html>
<!-- Created by images2html Erin Sheldon -->
<head>\n";
    <link rel="STYLESHEET" type="text/css" href="{css_url}>
</head>
    
<body>
    
    <div id="content">

        <div class="imagetable">
            <em>{image_name}</em><br>

	        <a href="{next_image_url}"><img src="../.thumbnails/{reduced_name}"></a>

        </div>

        <div class="prevNext">
            <div class="prevNextRow">
                <div class="prevtd">
                    {prev_image_url}
                </div>
                <div class="firstlast">
                    {index_url}
                </div>
                <div class="nexttd">
                    {next_image_url}
                </div>
            </div>
        </div>

    </div>

"""
