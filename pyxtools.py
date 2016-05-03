from __future__ import print_function
import os
import numpy
import tempfile
from pyx import *
from pyx.graph import axis


def plot(x, y, dx=None, dy=None, **kw):
    """
    first try to make a sensible interactive wrapper

    parameters
    ----------
    x: array or list
        x data to plot
    y: array or list
        y data to plot
    dx: array or list
        Error bars on x
    dy: array or list
        Error bars on y
    sym: string, optional
        Symbol name (circle,cross,diamond,plus,square,triangle)
    line: string, optional
        Line type
    color: string
        Color name from rgb.txt
    g: graphxy object
        Plot data into this object instead of creating a new one
    file: string
        File to write
    show: bool
        If True, bring up the image in a viewer

    keywords for graphxy
    """

    plt=kw.pop('plt',None)

    if plt is None:
        plt=Plotter(**kw)

    plt.plot(x, y, dx=dx, dy=dy, **kw)

    _writefile_maybe(plt.g, **kw)
    _show_maybe(plt.g, **kw)

    return plt

def imview(image, **kw):
    """
    plot the image

    parameters
    ----------
    image: array
        A numpy array
    """
    if len(image.shape) != 2:
        raise ValueError("image should be 2d")

    # we can allow the user to specify physical bounds as well
    xmin = 0
    xmax = image.shape[0]-1
    ymin = 0
    ymax = image.shape[1]-1
    x, y = numpy.mgrid[
        0:image.shape[0],
        0:image.shape[1],
    ]


    # need to convert to lists for pyx
    data = list(zip(x.flat, y.flat, image.flat))

    kw['xmin']=xmin
    kw['xmax']=xmax
    kw['ymin']=ymin
    kw['ymax']=ymax
    xaxis, xlog, yaxis, ylog=_get_axes(kw)
    assert xlog==False and ylog==False,"no log axes for images"

    g = graph.graphxy(
        height=8, width=8,
        x=xaxis,
        y=yaxis,
    )

    scale_title=kw.get('scale_title','')
    coloraxis=graph.axis.linear(
        min=0,
        max=image.max(),
        title=scale_title
    )

    # make a keyword for the color scheme
    gradient=color.gradient.ReverseGrey
    pstyle=graph.style.density(gradient=gradient, coloraxis=coloraxis)

    pdata=graph.data.points(data, x=1, y=2, color=3)
    g.plot(pdata, styles=[pstyle])

    _writefile_maybe(g, **kw)
    _show_maybe(g, **kw)

    return g


def write(g, fname, **kw):
    """
    Write the pyx object to a file. Supported types are
    png, jpg, pdf, eps, ps.  For newer pyx svg is supported

    parameters
    ----------
    pyxobj: a pyx object
        e.g. a graph or canvas
    fname: string
        File name.  The type will be inferred from the extensions
    **kw:
        other keywords such as resolution (or dpi for short)
    """
    _writefile_maybe(g, file=fname, **kw)

def _writefile_maybe(g, **kw):
    fname=kw.get('file',None)
    if fname is not None:
        if 'png' in fname or 'jpg' in fname:

            if 'dpi' in kw:
                res=kw['dpi']
            else:
                res=kw.get('resolution',200)

            g.writeGSfile(fname, resolution=res)
        else:
            g.writetofile(fname)

def _show_maybe(g, **kw):
    if True==kw.get('show',False):
        _do_show(g, **kw)


def _do_show(g, **kw):
    import time
    fname = tempfile.mktemp(suffix='.png')



    kw['file'] = fname
    _writefile_maybe(g, **kw)

    _tflist.add(fname)


    viewer=_config['viewer']
    #cmd='{viewer} {fname} &> /dev/null &'
    cmd='{viewer} {fname}  &'
    cmd=cmd.format(viewer=viewer,fname=fname)
    os.system(cmd)


class Plotter(object):
    """
    for use with the interactive plot() command
    """
    def __init__(self, **kw):
        self.value_sets=[]
        self.style_sets=[]

        self.graph_kw=kw
        self.xlog=kw.get('xlog',False)
        self.ylog=kw.get('ylog',False)
        self.g=None

        self.xmin=None
        self.xmax=None
        self.ymin=None
        self.ymax=None

    def plot(self, x, y, dx=None, dy=None, **kw):
        """
        add values to the plot
        """
        styles=[]

        values, xrng, yrng =self._get_values(x, y, dx=dx, dy=dy)


        self._set_ranges(xrng, yrng)

        #if self.g is None:
        #    self.g=self._get_graph_and_axes(xrng, yrng, **self.graph_kw)
        self.g=self._get_graph_and_axes(xrng, yrng, **self.graph_kw)

        styles=self._set_symbol(styles, dx=dx, dy=dy, **kw)

        self.value_sets.append(values)
        self.style_sets.append(styles)

        for i in xrange(len(self.value_sets)):
            self.g.plot(
                self.value_sets[i],
                styles=self.style_sets[i]
            )

    def write(self, filename, **kw):
        """
        write to a plot file

        parameters
        ----------
        filename: string
            file type determined from extension
        **kw:
            extra keywords
        """
        if self.g is None:
            raise RuntimeError("plot some data first")

        kw['file'] = filename
        _writefile_maybe(self.g, **kw)

    def show(self, **kw):
        """
        show the plot

        parameters
        ----------
        **kw:
            extra keywords
        """
        if self.g is None:
            raise RuntimeError("plot some data first")

        kw['show']=True
        _show_maybe(self.g, **kw)

    def _set_ranges(self, xrng, yrng):
        if self.xmin is None:
            self.xmin=xrng[0]
        else:
            self.xmin=min(self.xmin, xrng[0])

        if self.xmax is None:
            self.xmax=xrng[1]
        else:
            self.xmax=max(self.xmax, xrng[1])


        if self.ymin is None:
            self.ymin=yrng[0]
        else:
            self.ymin=min(self.ymin, yrng[0])

        if self.ymax is None:
            self.ymax=yrng[1]
        else:
            self.ymax=max(self.ymax, yrng[1])


    def _get_graph_and_axes(self, xrng, yrng, **kw):
        gkw=_unpack_graphxy_keywords(kw)
        gkw['width'] = gkw.get('width',8)

        self.xlog=kw.get('xlog',False)
        self.ylog=kw.get('ylog',False)

        xdiff=xrng[1]-xrng[0]
        ydiff=yrng[1]-yrng[0]

        if 'xmin' not in kw:
            kw['xmin'] = _get_prng(xrng[0], xdiff, 'low', log=self.xlog)
        if 'xmax' not in kw:
            kw['xmax'] = _get_prng(xrng[1], xdiff, 'high', log=self.xlog)
        if 'ymin' not in kw:
            kw['ymin'] = _get_prng(yrng[0], ydiff, 'low', log=self.ylog)
        if 'ymax' not in kw:
            kw['ymax'] = _get_prng(yrng[1], ydiff, 'high', log=self.ylog)

        xaxis,xlog,yaxis,ylog = _get_axes(kw)
        g = graph.graphxy(x=xaxis, y=yaxis, **gkw)

        return g

    def _set_symbol(self, styles, dx=None, dy=None, **kw):

        if 'sym' in kw:
            sym=symbols.get_symbol(kw['sym'])
        else:
            sym=symbols.get_symbol('circle')

        if 'color' in kw:
            clr=colors(kw['color'])
        else:
            clr=colors('black')

        symbolattrs=[clr,deco.filled([clr])]

        symbol=graph.style.symbol(
            symbol=sym,
            size=kw.get('size',0.1),
            symbolattrs=symbolattrs,
        )
        styles += [symbol]

        if dx is not None or dy is not None:
            if 'errcolor' in kw:
                errclr=colors(kw['errcolor'])
            else:
                errclr=clr
            styles += [graph.style.errorbar(errorbarattrs=[errclr])]

        return styles

    def _get_values(self, x, y, dx=None, dy=None):
        args={'x':list(x), 'y':list(y)}

        xrng=[x.min(), x.max()]
        yrng=[y.min(), y.max()]

        if dx is not None:
            if self.xlog:
                # clip way below the value
                w,=numpy.where(x > 0)
                if w.size > 0:
                    # need something better than this
                    lowest=1.0e-6*x[w].min()
                    xlower=(x-dx).clip(lowest)
                    xupper=(x+dx).clip(lowest)

                    args['xmin'] = xlower
                    args['xmax'] = xupper

                    xrng[0]=xlower.min()
                    xrng[1]=xupper.max()
            else:
                xrng[0] = (x-dx).min()
                xrng[1] = (x+dx).max()
                args['dx'] = dx

        if dy is not None:
            if self.ylog:
                # clip way below the value
                w,=numpy.where(y > 0)
                if w.size > 0:
                    # need something better than this
                    lowest=1.0e-6*y[w].min()
                    ylower=(y-dy).clip(lowest)
                    yupper=(y+dy).clip(lowest)

                    args['ymin'] = ylower
                    args['ymax'] = yupper

                    yrng[0]=ylower.min()
                    yrng[1]=yupper.max()
            else:
                yrng[0] = (y-dy).min()
                yrng[1] = (y+dy).max()
                args['dy'] = dy


        values=graph.data.values(**args)

        return values, xrng, yrng

_graphxy_kw=[
    'xpos','ypos','width', 'height', 'ratio',
    'key','backgroundattrs','axesdist', 'flipped',
    'xaxisat', 'yaxisat', 'axes']

def _get_prngold(x, type, log=False):
    if log:
        if type=='low':
            r=0.5*x
        else:
            r=1.5*x
    else:
        if type=='low':
            if x < 0:
                r=1.2*x
            else:
                r=0.8*x
        else:
            if x < 0:
                r=0.8*x
            else:
                r=1.2*x

    return r

def _get_prng(x, diff, type, log=False):

    frac=0.075
    if log:
        if type=='low':
            r=0.5*x
        else:
            r=1.5*x
    else:

        fdiff=diff*frac
        if type=='low':
            r = x - fdiff
        else:
            r = x + fdiff

    print("x:",x,"r:",r)
    return r


def _unpack_graphxy_keywords(kwin):
    kw={}
    for key in kwin:
        if key in _graphxy_kw:
            kw[key] = kwin[key]

    return kw

def _get_axes(kw):
    xlog=kw.get('xlog',False)
    ylog=kw.get('ylog',False)

    xkw=dict(
        title=kw.get('xlabel',None),

        min=kw.get('xmin',None),
        max=kw.get('xmax',None),

        # density of ticks
        density=kw.get('xdensity',1),
        reverse=kw.get('xreverse',0),
    )
    ykw=dict(
        title=kw.get('ylabel',None),

        min=kw.get('ymin',None),
        max=kw.get('ymax',None),

        # density of ticks
        density=kw.get('ydensity',1),
        reverse=kw.get('yreverse',0),
    )

    if xlog:
        xaxis = axis.log(**xkw)
    else:
        xaxis = axis.lin(**xkw)

    if ylog:
        yaxis = axis.log(**ykw)
    else:
        yaxis = axis.lin(**ykw)

    return xaxis, xlog, yaxis, ylog

def test_log():
    import numpy
    from numpy import log10
    x=numpy.logspace(log10(0.1), log10(30.0), 10)
    y=1.0/x

    # we here use parters and texters which are explained in the examples below
    #log2parter = axis.parter.log([axis.parter.preexp([axis.tick.rational(1)], 4),
    #                              axis.parter.preexp([axis.tick.rational(1)], 2)])
    #log2texter = axis.texter.exponential(nomantissaexp=r"{2^{%s}}",
    #                                     mantissamax=axis.tick.rational(2))

    fac=2.0
    g = graph.graphxy(
        width=8,
        height=8,
        #ratio=1.0,
        x=axis.log(min=x.min()/fac,max=fac*x.max(),title=r'$R~ [h^{-1}$ Mpc]'),
        y=axis.log(min=y.min()/fac,max=fac*y.max(),title=r'$\Delta\Sigma ~[h~$M$_\odot~$pc$^{-2}]$'),
    )

    symbol=graph.style.symbol(
        symbol=graph.style._circlesymbol,
        size=0.1,
        symbolattrs=[deco.filled([color.rgb.black])],
    )

    values=graph.data.values(x=list(x),y=list(y))
    g.plot(values,[symbol])

    #g.writeEPSfile("log")
    #g.writePDFfile("log")
    g.writetofile("log-test.pdf")
    g.writetofile("log-test.eps")

def test_basic():

    # can define by hex string or rgb value, normalized to 1
    #yellow=color.rgbfromhexstring('#FFFF00')
    yellow=color.rgb(r=1.0, g=1.0, b=0.0)

    # need to specify one or both of width, height
    #key=graph.key.key(pos="br", dist=0.1)
    key=graph.key.key(pos='tl')
    g = graph.graphxy(
        width=8,
        height=8,
        key=key,
        x=axis.lin(title=r"$x$-axis"),
        y=axis.lin(title=r"$y$-axis"),
    )


    # either provide lists of the individual coordinates
    factor=1.2343
    x=list(range(10))
    y=[factor*xx**2 for xx in x]
    values=graph.data.values(x=x, y=y,title='points')
    symbol1=graph.style.symbol(
        #symbol=graph.style._circlesymbol,
        symbol=graph.style.symbol.circle,
        symbolattrs=[deco.filled([yellow])],
    )
    symbol2=graph.style.symbol(
        #symbol=graph.style._circlesymbol,
        symbol=graph.style.symbol.circle,
        symbolattrs=[color.rgb.red],
    )

    g.plot(values,[symbol1,symbol2])

    # or provide one list containing the whole points
    #values2=graph.data.points(list(zip(range(10), range(10))), x=1, y=2)
    #g.plot(values2)

    # also a function
    c=graph.data.function(
        "y(x)=%g * x**2" % factor,
        min=0, max=9,
        title=r'$y=x^2$'
    )
    g.plot(c)

    g.writeEPSfile("points-and-line")
    g.writePDFfile("points-and-line")
    #g.writeSVGfile("points")

def test_image(**kw):
    import numpy
    from pyx import color, graph

    xmax = 1.6
    xmin = -xmax
    ymax = 1.6
    ymin = -ymax
    npts = 101
    sigma=0.5
    x, y = numpy.mgrid[xmin:xmax:npts*1j, ymin:ymax:npts*1j]
    z = numpy.exp( -(x**2 + y**2)/(2*sigma**2) )

    g=imview(z, **kw)
    return g

_symdict={
    'circle': graph.style.symbol.circle,
    'cross': graph.style.symbol.cross,
    'diamond': graph.style.symbol.diamond,
    'plus': graph.style.symbol.plus,
    'square': graph.style.symbol.square,
    'triangle': graph.style.symbol.triangle,
}

class Symbols(object):
    def get_symbol(self, symname):
        if isinstance(symname,  graph.style.symbol):
            return symname

        sym=_symdict.get(symname,None)
        if sym is None:
            raise ValueError("bad symbol name: '%s'" % symname)
        return sym

symbols=Symbols()

_linedict={
    'solid': style.linestyle.solid,
}

class Linestyles(object):
    def get_style(self, linename):
        line=_linedict.get(linename,None)
        if line is None:
            raise ValueError("bad line style name: '%s'" % linename)
        return line

linestyles=Linestyles()


class Color(object):
    """
    translate color names from rgb.txt to pyx.color.rgb instances
    """
    def __call__(self, name):
        return self.get_rgb(name)

    def get_rgb(self, name_in):
        if isinstance(name_in, color.rgb):
            return name_in

        name=name_in.lower()
        rgb=_rgbdict.get(name,None)
        if rgb is None:
            raise ValueError("bad color name: '%s'" % name_in)

        r=rgb[0]/255.0
        g=rgb[1]/255.0
        b=rgb[2]/255.0
        return color.rgb(r=r, g=g, b=b)

colors=Color()



def _make_rgb_colors_dict(fname):
    """
    fname is the rgb.txt file

    255 228 196             bisque
    """

    rgbdict={}
    with open(fname) as fobj:
        for i,line in enumerate(fobj):
            if i==0:
                continue

            ls=line.split()

            # skip names with blank spaces
            if len(ls) > 4:
                continue

            r=int(ls[0])
            g=int(ls[1])
            b=int(ls[2])
            name=ls[3].lower()

            rgbdict[name] = (r,g,b)

    return rgbdict


_rgbdict = {
    'aliceblue': (240, 248, 255),
    'antiquewhite': (250, 235, 215),
    'antiquewhite1': (255, 239, 219),
    'antiquewhite2': (238, 223, 204),
    'antiquewhite3': (205, 192, 176),
    'antiquewhite4': (139, 131, 120),
    'aquamarine': (127, 255, 212),
    'aquamarine1': (127, 255, 212),
    'aquamarine2': (118, 238, 198),
    'aquamarine3': (102, 205, 170),
    'aquamarine4': (69, 139, 116),
    'azure': (240, 255, 255),
    'azure1': (240, 255, 255),
    'azure2': (224, 238, 238),
    'azure3': (193, 205, 205),
    'azure4': (131, 139, 139),
    'beige': (245, 245, 220),
    'bisque': (255, 228, 196),
    'bisque1': (255, 228, 196),
    'bisque2': (238, 213, 183),
    'bisque3': (205, 183, 158),
    'bisque4': (139, 125, 107),
    'black': (0, 0, 0),
    'blanchedalmond': (255, 235, 205),
    'blue': (0, 0, 255),
    'blue1': (0, 0, 255),
    'blue2': (0, 0, 238),
    'blue3': (0, 0, 205),
    'blue4': (0, 0, 139),
    'blueviolet': (138, 43, 226),
    'brown': (165, 42, 42),
    'brown1': (255, 64, 64),
    'brown2': (238, 59, 59),
    'brown3': (205, 51, 51),
    'brown4': (139, 35, 35),
    'burlywood': (222, 184, 135),
    'burlywood1': (255, 211, 155),
    'burlywood2': (238, 197, 145),
    'burlywood3': (205, 170, 125),
    'burlywood4': (139, 115, 85),
    'cadetblue': (95, 158, 160),
    'cadetblue1': (152, 245, 255),
    'cadetblue2': (142, 229, 238),
    'cadetblue3': (122, 197, 205),
    'cadetblue4': (83, 134, 139),
    'chartreuse': (127, 255, 0),
    'chartreuse1': (127, 255, 0),
    'chartreuse2': (118, 238, 0),
    'chartreuse3': (102, 205, 0),
    'chartreuse4': (69, 139, 0),
    'chocolate': (210, 105, 30),
    'chocolate1': (255, 127, 36),
    'chocolate2': (238, 118, 33),
    'chocolate3': (205, 102, 29),
    'chocolate4': (139, 69, 19),
    'coral': (255, 127, 80),
    'coral1': (255, 114, 86),
    'coral2': (238, 106, 80),
    'coral3': (205, 91, 69),
    'coral4': (139, 62, 47),
    'cornflowerblue': (100, 149, 237),
    'cornsilk': (255, 248, 220),
    'cornsilk1': (255, 248, 220),
    'cornsilk2': (238, 232, 205),
    'cornsilk3': (205, 200, 177),
    'cornsilk4': (139, 136, 120),
    'cyan': (0, 255, 255),
    'cyan1': (0, 255, 255),
    'cyan2': (0, 238, 238),
    'cyan3': (0, 205, 205),
    'cyan4': (0, 139, 139),
    'darkblue': (0, 0, 139),
    'darkcyan': (0, 139, 139),
    'darkgoldenrod': (184, 134, 11),
    'darkgoldenrod1': (255, 185, 15),
    'darkgoldenrod2': (238, 173, 14),
    'darkgoldenrod3': (205, 149, 12),
    'darkgoldenrod4': (139, 101, 8),
    'darkgray': (169, 169, 169),
    'darkgreen': (0, 100, 0),
    'darkgrey': (169, 169, 169),
    'darkkhaki': (189, 183, 107),
    'darkmagenta': (139, 0, 139),
    'darkolivegreen': (85, 107, 47),
    'darkolivegreen1': (202, 255, 112),
    'darkolivegreen2': (188, 238, 104),
    'darkolivegreen3': (162, 205, 90),
    'darkolivegreen4': (110, 139, 61),
    'darkorange': (255, 140, 0),
    'darkorange1': (255, 127, 0),
    'darkorange2': (238, 118, 0),
    'darkorange3': (205, 102, 0),
    'darkorange4': (139, 69, 0),
    'darkorchid': (153, 50, 204),
    'darkorchid1': (191, 62, 255),
    'darkorchid2': (178, 58, 238),
    'darkorchid3': (154, 50, 205),
    'darkorchid4': (104, 34, 139),
    'darkred': (139, 0, 0),
    'darksalmon': (233, 150, 122),
    'darkseagreen': (143, 188, 143),
    'darkseagreen1': (193, 255, 193),
    'darkseagreen2': (180, 238, 180),
    'darkseagreen3': (155, 205, 155),
    'darkseagreen4': (105, 139, 105),
    'darkslateblue': (72, 61, 139),
    'darkslategray': (47, 79, 79),
    'darkslategray1': (151, 255, 255),
    'darkslategray2': (141, 238, 238),
    'darkslategray3': (121, 205, 205),
    'darkslategray4': (82, 139, 139),
    'darkslategrey': (47, 79, 79),
    'darkturquoise': (0, 206, 209),
    'darkviolet': (148, 0, 211),
    'debianred': (215, 7, 81),
    'deeppink': (255, 20, 147),
    'deeppink1': (255, 20, 147),
    'deeppink2': (238, 18, 137),
    'deeppink3': (205, 16, 118),
    'deeppink4': (139, 10, 80),
    'deepskyblue': (0, 191, 255),
    'deepskyblue1': (0, 191, 255),
    'deepskyblue2': (0, 178, 238),
    'deepskyblue3': (0, 154, 205),
    'deepskyblue4': (0, 104, 139),
    'dimgray': (105, 105, 105),
    'dimgrey': (105, 105, 105),
    'dodgerblue': (30, 144, 255),
    'dodgerblue1': (30, 144, 255),
    'dodgerblue2': (28, 134, 238),
    'dodgerblue3': (24, 116, 205),
    'dodgerblue4': (16, 78, 139),
    'firebrick': (178, 34, 34),
    'firebrick1': (255, 48, 48),
    'firebrick2': (238, 44, 44),
    'firebrick3': (205, 38, 38),
    'firebrick4': (139, 26, 26),
    'floralwhite': (255, 250, 240),
    'forestgreen': (34, 139, 34),
    'gainsboro': (220, 220, 220),
    'ghostwhite': (248, 248, 255),
    'gold': (255, 215, 0),
    'gold1': (255, 215, 0),
    'gold2': (238, 201, 0),
    'gold3': (205, 173, 0),
    'gold4': (139, 117, 0),
    'goldenrod': (218, 165, 32),
    'goldenrod1': (255, 193, 37),
    'goldenrod2': (238, 180, 34),
    'goldenrod3': (205, 155, 29),
    'goldenrod4': (139, 105, 20),
    'gray': (190, 190, 190),
    'gray0': (0, 0, 0),
    'gray1': (3, 3, 3),
    'gray10': (26, 26, 26),
    'gray100': (255, 255, 255),
    'gray11': (28, 28, 28),
    'gray12': (31, 31, 31),
    'gray13': (33, 33, 33),
    'gray14': (36, 36, 36),
    'gray15': (38, 38, 38),
    'gray16': (41, 41, 41),
    'gray17': (43, 43, 43),
    'gray18': (46, 46, 46),
    'gray19': (48, 48, 48),
    'gray2': (5, 5, 5),
    'gray20': (51, 51, 51),
    'gray21': (54, 54, 54),
    'gray22': (56, 56, 56),
    'gray23': (59, 59, 59),
    'gray24': (61, 61, 61),
    'gray25': (64, 64, 64),
    'gray26': (66, 66, 66),
    'gray27': (69, 69, 69),
    'gray28': (71, 71, 71),
    'gray29': (74, 74, 74),
    'gray3': (8, 8, 8),
    'gray30': (77, 77, 77),
    'gray31': (79, 79, 79),
    'gray32': (82, 82, 82),
    'gray33': (84, 84, 84),
    'gray34': (87, 87, 87),
    'gray35': (89, 89, 89),
    'gray36': (92, 92, 92),
    'gray37': (94, 94, 94),
    'gray38': (97, 97, 97),
    'gray39': (99, 99, 99),
    'gray4': (10, 10, 10),
    'gray40': (102, 102, 102),
    'gray41': (105, 105, 105),
    'gray42': (107, 107, 107),
    'gray43': (110, 110, 110),
    'gray44': (112, 112, 112),
    'gray45': (115, 115, 115),
    'gray46': (117, 117, 117),
    'gray47': (120, 120, 120),
    'gray48': (122, 122, 122),
    'gray49': (125, 125, 125),
    'gray5': (13, 13, 13),
    'gray50': (127, 127, 127),
    'gray51': (130, 130, 130),
    'gray52': (133, 133, 133),
    'gray53': (135, 135, 135),
    'gray54': (138, 138, 138),
    'gray55': (140, 140, 140),
    'gray56': (143, 143, 143),
    'gray57': (145, 145, 145),
    'gray58': (148, 148, 148),
    'gray59': (150, 150, 150),
    'gray6': (15, 15, 15),
    'gray60': (153, 153, 153),
    'gray61': (156, 156, 156),
    'gray62': (158, 158, 158),
    'gray63': (161, 161, 161),
    'gray64': (163, 163, 163),
    'gray65': (166, 166, 166),
    'gray66': (168, 168, 168),
    'gray67': (171, 171, 171),
    'gray68': (173, 173, 173),
    'gray69': (176, 176, 176),
    'gray7': (18, 18, 18),
    'gray70': (179, 179, 179),
    'gray71': (181, 181, 181),
    'gray72': (184, 184, 184),
    'gray73': (186, 186, 186),
    'gray74': (189, 189, 189),
    'gray75': (191, 191, 191),
    'gray76': (194, 194, 194),
    'gray77': (196, 196, 196),
    'gray78': (199, 199, 199),
    'gray79': (201, 201, 201),
    'gray8': (20, 20, 20),
    'gray80': (204, 204, 204),
    'gray81': (207, 207, 207),
    'gray82': (209, 209, 209),
    'gray83': (212, 212, 212),
    'gray84': (214, 214, 214),
    'gray85': (217, 217, 217),
    'gray86': (219, 219, 219),
    'gray87': (222, 222, 222),
    'gray88': (224, 224, 224),
    'gray89': (227, 227, 227),
    'gray9': (23, 23, 23),
    'gray90': (229, 229, 229),
    'gray91': (232, 232, 232),
    'gray92': (235, 235, 235),
    'gray93': (237, 237, 237),
    'gray94': (240, 240, 240),
    'gray95': (242, 242, 242),
    'gray96': (245, 245, 245),
    'gray97': (247, 247, 247),
    'gray98': (250, 250, 250),
    'gray99': (252, 252, 252),
    'green': (0, 255, 0),
    'green1': (0, 255, 0),
    'green2': (0, 238, 0),
    'green3': (0, 205, 0),
    'green4': (0, 139, 0),
    'greenyellow': (173, 255, 47),
    'grey': (190, 190, 190),
    'grey0': (0, 0, 0),
    'grey1': (3, 3, 3),
    'grey10': (26, 26, 26),
    'grey100': (255, 255, 255),
    'grey11': (28, 28, 28),
    'grey12': (31, 31, 31),
    'grey13': (33, 33, 33),
    'grey14': (36, 36, 36),
    'grey15': (38, 38, 38),
    'grey16': (41, 41, 41),
    'grey17': (43, 43, 43),
    'grey18': (46, 46, 46),
    'grey19': (48, 48, 48),
    'grey2': (5, 5, 5),
    'grey20': (51, 51, 51),
    'grey21': (54, 54, 54),
    'grey22': (56, 56, 56),
    'grey23': (59, 59, 59),
    'grey24': (61, 61, 61),
    'grey25': (64, 64, 64),
    'grey26': (66, 66, 66),
    'grey27': (69, 69, 69),
    'grey28': (71, 71, 71),
    'grey29': (74, 74, 74),
    'grey3': (8, 8, 8),
    'grey30': (77, 77, 77),
    'grey31': (79, 79, 79),
    'grey32': (82, 82, 82),
    'grey33': (84, 84, 84),
    'grey34': (87, 87, 87),
    'grey35': (89, 89, 89),
    'grey36': (92, 92, 92),
    'grey37': (94, 94, 94),
    'grey38': (97, 97, 97),
    'grey39': (99, 99, 99),
    'grey4': (10, 10, 10),
    'grey40': (102, 102, 102),
    'grey41': (105, 105, 105),
    'grey42': (107, 107, 107),
    'grey43': (110, 110, 110),
    'grey44': (112, 112, 112),
    'grey45': (115, 115, 115),
    'grey46': (117, 117, 117),
    'grey47': (120, 120, 120),
    'grey48': (122, 122, 122),
    'grey49': (125, 125, 125),
    'grey5': (13, 13, 13),
    'grey50': (127, 127, 127),
    'grey51': (130, 130, 130),
    'grey52': (133, 133, 133),
    'grey53': (135, 135, 135),
    'grey54': (138, 138, 138),
    'grey55': (140, 140, 140),
    'grey56': (143, 143, 143),
    'grey57': (145, 145, 145),
    'grey58': (148, 148, 148),
    'grey59': (150, 150, 150),
    'grey6': (15, 15, 15),
    'grey60': (153, 153, 153),
    'grey61': (156, 156, 156),
    'grey62': (158, 158, 158),
    'grey63': (161, 161, 161),
    'grey64': (163, 163, 163),
    'grey65': (166, 166, 166),
    'grey66': (168, 168, 168),
    'grey67': (171, 171, 171),
    'grey68': (173, 173, 173),
    'grey69': (176, 176, 176),
    'grey7': (18, 18, 18),
    'grey70': (179, 179, 179),
    'grey71': (181, 181, 181),
    'grey72': (184, 184, 184),
    'grey73': (186, 186, 186),
    'grey74': (189, 189, 189),
    'grey75': (191, 191, 191),
    'grey76': (194, 194, 194),
    'grey77': (196, 196, 196),
    'grey78': (199, 199, 199),
    'grey79': (201, 201, 201),
    'grey8': (20, 20, 20),
    'grey80': (204, 204, 204),
    'grey81': (207, 207, 207),
    'grey82': (209, 209, 209),
    'grey83': (212, 212, 212),
    'grey84': (214, 214, 214),
    'grey85': (217, 217, 217),
    'grey86': (219, 219, 219),
    'grey87': (222, 222, 222),
    'grey88': (224, 224, 224),
    'grey89': (227, 227, 227),
    'grey9': (23, 23, 23),
    'grey90': (229, 229, 229),
    'grey91': (232, 232, 232),
    'grey92': (235, 235, 235),
    'grey93': (237, 237, 237),
    'grey94': (240, 240, 240),
    'grey95': (242, 242, 242),
    'grey96': (245, 245, 245),
    'grey97': (247, 247, 247),
    'grey98': (250, 250, 250),
    'grey99': (252, 252, 252),
    'honeydew': (240, 255, 240),
    'honeydew1': (240, 255, 240),
    'honeydew2': (224, 238, 224),
    'honeydew3': (193, 205, 193),
    'honeydew4': (131, 139, 131),
    'hotpink': (255, 105, 180),
    'hotpink1': (255, 110, 180),
    'hotpink2': (238, 106, 167),
    'hotpink3': (205, 96, 144),
    'hotpink4': (139, 58, 98),
    'indianred': (205, 92, 92),
    'indianred1': (255, 106, 106),
    'indianred2': (238, 99, 99),
    'indianred3': (205, 85, 85),
    'indianred4': (139, 58, 58),
    'ivory': (255, 255, 240),
    'ivory1': (255, 255, 240),
    'ivory2': (238, 238, 224),
    'ivory3': (205, 205, 193),
    'ivory4': (139, 139, 131),
    'khaki': (240, 230, 140),
    'khaki1': (255, 246, 143),
    'khaki2': (238, 230, 133),
    'khaki3': (205, 198, 115),
    'khaki4': (139, 134, 78),
    'lavender': (230, 230, 250),
    'lavenderblush': (255, 240, 245),
    'lavenderblush1': (255, 240, 245),
    'lavenderblush2': (238, 224, 229),
    'lavenderblush3': (205, 193, 197),
    'lavenderblush4': (139, 131, 134),
    'lawngreen': (124, 252, 0),
    'lemonchiffon': (255, 250, 205),
    'lemonchiffon1': (255, 250, 205),
    'lemonchiffon2': (238, 233, 191),
    'lemonchiffon3': (205, 201, 165),
    'lemonchiffon4': (139, 137, 112),
    'lightblue': (173, 216, 230),
    'lightblue1': (191, 239, 255),
    'lightblue2': (178, 223, 238),
    'lightblue3': (154, 192, 205),
    'lightblue4': (104, 131, 139),
    'lightcoral': (240, 128, 128),
    'lightcyan': (224, 255, 255),
    'lightcyan1': (224, 255, 255),
    'lightcyan2': (209, 238, 238),
    'lightcyan3': (180, 205, 205),
    'lightcyan4': (122, 139, 139),
    'lightgoldenrod': (238, 221, 130),
    'lightgoldenrod1': (255, 236, 139),
    'lightgoldenrod2': (238, 220, 130),
    'lightgoldenrod3': (205, 190, 112),
    'lightgoldenrod4': (139, 129, 76),
    'lightgoldenrodyellow': (250, 250, 210),
    'lightgray': (211, 211, 211),
    'lightgreen': (144, 238, 144),
    'lightgrey': (211, 211, 211),
    'lightpink': (255, 182, 193),
    'lightpink1': (255, 174, 185),
    'lightpink2': (238, 162, 173),
    'lightpink3': (205, 140, 149),
    'lightpink4': (139, 95, 101),
    'lightsalmon': (255, 160, 122),
    'lightsalmon1': (255, 160, 122),
    'lightsalmon2': (238, 149, 114),
    'lightsalmon3': (205, 129, 98),
    'lightsalmon4': (139, 87, 66),
    'lightseagreen': (32, 178, 170),
    'lightskyblue': (135, 206, 250),
    'lightskyblue1': (176, 226, 255),
    'lightskyblue2': (164, 211, 238),
    'lightskyblue3': (141, 182, 205),
    'lightskyblue4': (96, 123, 139),
    'lightslateblue': (132, 112, 255),
    'lightslategray': (119, 136, 153),
    'lightslategrey': (119, 136, 153),
    'lightsteelblue': (176, 196, 222),
    'lightsteelblue1': (202, 225, 255),
    'lightsteelblue2': (188, 210, 238),
    'lightsteelblue3': (162, 181, 205),
    'lightsteelblue4': (110, 123, 139),
    'lightyellow': (255, 255, 224),
    'lightyellow1': (255, 255, 224),
    'lightyellow2': (238, 238, 209),
    'lightyellow3': (205, 205, 180),
    'lightyellow4': (139, 139, 122),
    'limegreen': (50, 205, 50),
    'linen': (250, 240, 230),
    'magenta': (255, 0, 255),
    'magenta1': (255, 0, 255),
    'magenta2': (238, 0, 238),
    'magenta3': (205, 0, 205),
    'magenta4': (139, 0, 139),
    'maroon': (176, 48, 96),
    'maroon1': (255, 52, 179),
    'maroon2': (238, 48, 167),
    'maroon3': (205, 41, 144),
    'maroon4': (139, 28, 98),
    'mediumaquamarine': (102, 205, 170),
    'mediumblue': (0, 0, 205),
    'mediumorchid': (186, 85, 211),
    'mediumorchid1': (224, 102, 255),
    'mediumorchid2': (209, 95, 238),
    'mediumorchid3': (180, 82, 205),
    'mediumorchid4': (122, 55, 139),
    'mediumpurple': (147, 112, 219),
    'mediumpurple1': (171, 130, 255),
    'mediumpurple2': (159, 121, 238),
    'mediumpurple3': (137, 104, 205),
    'mediumpurple4': (93, 71, 139),
    'mediumseagreen': (60, 179, 113),
    'mediumslateblue': (123, 104, 238),
    'mediumspringgreen': (0, 250, 154),
    'mediumturquoise': (72, 209, 204),
    'mediumvioletred': (199, 21, 133),
    'midnightblue': (25, 25, 112),
    'mintcream': (245, 255, 250),
    'mistyrose': (255, 228, 225),
    'mistyrose1': (255, 228, 225),
    'mistyrose2': (238, 213, 210),
    'mistyrose3': (205, 183, 181),
    'mistyrose4': (139, 125, 123),
    'moccasin': (255, 228, 181),
    'navajowhite': (255, 222, 173),
    'navajowhite1': (255, 222, 173),
    'navajowhite2': (238, 207, 161),
    'navajowhite3': (205, 179, 139),
    'navajowhite4': (139, 121, 94),
    'navy': (0, 0, 128),
    'navyblue': (0, 0, 128),
    'oldlace': (253, 245, 230),
    'olivedrab': (107, 142, 35),
    'olivedrab1': (192, 255, 62),
    'olivedrab2': (179, 238, 58),
    'olivedrab3': (154, 205, 50),
    'olivedrab4': (105, 139, 34),
    'orange': (255, 165, 0),
    'orange1': (255, 165, 0),
    'orange2': (238, 154, 0),
    'orange3': (205, 133, 0),
    'orange4': (139, 90, 0),
    'orangered': (255, 69, 0),
    'orangered1': (255, 69, 0),
    'orangered2': (238, 64, 0),
    'orangered3': (205, 55, 0),
    'orangered4': (139, 37, 0),
    'orchid': (218, 112, 214),
    'orchid1': (255, 131, 250),
    'orchid2': (238, 122, 233),
    'orchid3': (205, 105, 201),
    'orchid4': (139, 71, 137),
    'palegoldenrod': (238, 232, 170),
    'palegreen': (152, 251, 152),
    'palegreen1': (154, 255, 154),
    'palegreen2': (144, 238, 144),
    'palegreen3': (124, 205, 124),
    'palegreen4': (84, 139, 84),
    'paleturquoise': (175, 238, 238),
    'paleturquoise1': (187, 255, 255),
    'paleturquoise2': (174, 238, 238),
    'paleturquoise3': (150, 205, 205),
    'paleturquoise4': (102, 139, 139),
    'palevioletred': (219, 112, 147),
    'palevioletred1': (255, 130, 171),
    'palevioletred2': (238, 121, 159),
    'palevioletred3': (205, 104, 137),
    'palevioletred4': (139, 71, 93),
    'papayawhip': (255, 239, 213),
    'peachpuff': (255, 218, 185),
    'peachpuff1': (255, 218, 185),
    'peachpuff2': (238, 203, 173),
    'peachpuff3': (205, 175, 149),
    'peachpuff4': (139, 119, 101),
    'peru': (205, 133, 63),
    'pink': (255, 192, 203),
    'pink1': (255, 181, 197),
    'pink2': (238, 169, 184),
    'pink3': (205, 145, 158),
    'pink4': (139, 99, 108),
    'plum': (221, 160, 221),
    'plum1': (255, 187, 255),
    'plum2': (238, 174, 238),
    'plum3': (205, 150, 205),
    'plum4': (139, 102, 139),
    'powderblue': (176, 224, 230),
    'purple': (160, 32, 240),
    'purple1': (155, 48, 255),
    'purple2': (145, 44, 238),
    'purple3': (125, 38, 205),
    'purple4': (85, 26, 139),
    'red': (255, 0, 0),
    'red1': (255, 0, 0),
    'red2': (238, 0, 0),
    'red3': (205, 0, 0),
    'red4': (139, 0, 0),
    'rosybrown': (188, 143, 143),
    'rosybrown1': (255, 193, 193),
    'rosybrown2': (238, 180, 180),
    'rosybrown3': (205, 155, 155),
    'rosybrown4': (139, 105, 105),
    'royalblue': (65, 105, 225),
    'royalblue1': (72, 118, 255),
    'royalblue2': (67, 110, 238),
    'royalblue3': (58, 95, 205),
    'royalblue4': (39, 64, 139),
    'saddlebrown': (139, 69, 19),
    'salmon': (250, 128, 114),
    'salmon1': (255, 140, 105),
    'salmon2': (238, 130, 98),
    'salmon3': (205, 112, 84),
    'salmon4': (139, 76, 57),
    'sandybrown': (244, 164, 96),
    'seagreen': (46, 139, 87),
    'seagreen1': (84, 255, 159),
    'seagreen2': (78, 238, 148),
    'seagreen3': (67, 205, 128),
    'seagreen4': (46, 139, 87),
    'seashell': (255, 245, 238),
    'seashell1': (255, 245, 238),
    'seashell2': (238, 229, 222),
    'seashell3': (205, 197, 191),
    'seashell4': (139, 134, 130),
    'sienna': (160, 82, 45),
    'sienna1': (255, 130, 71),
    'sienna2': (238, 121, 66),
    'sienna3': (205, 104, 57),
    'sienna4': (139, 71, 38),
    'skyblue': (135, 206, 235),
    'skyblue1': (135, 206, 255),
    'skyblue2': (126, 192, 238),
    'skyblue3': (108, 166, 205),
    'skyblue4': (74, 112, 139),
    'slateblue': (106, 90, 205),
    'slateblue1': (131, 111, 255),
    'slateblue2': (122, 103, 238),
    'slateblue3': (105, 89, 205),
    'slateblue4': (71, 60, 139),
    'slategray': (112, 128, 144),
    'slategray1': (198, 226, 255),
    'slategray2': (185, 211, 238),
    'slategray3': (159, 182, 205),
    'slategray4': (108, 123, 139),
    'slategrey': (112, 128, 144),
    'snow': (255, 250, 250),
    'snow1': (255, 250, 250),
    'snow2': (238, 233, 233),
    'snow3': (205, 201, 201),
    'snow4': (139, 137, 137),
    'springgreen': (0, 255, 127),
    'springgreen1': (0, 255, 127),
    'springgreen2': (0, 238, 118),
    'springgreen3': (0, 205, 102),
    'springgreen4': (0, 139, 69),
    'steelblue': (70, 130, 180),
    'steelblue1': (99, 184, 255),
    'steelblue2': (92, 172, 238),
    'steelblue3': (79, 148, 205),
    'steelblue4': (54, 100, 139),
    'tan': (210, 180, 140),
    'tan1': (255, 165, 79),
    'tan2': (238, 154, 73),
    'tan3': (205, 133, 63),
    'tan4': (139, 90, 43),
    'thistle': (216, 191, 216),
    'thistle1': (255, 225, 255),
    'thistle2': (238, 210, 238),
    'thistle3': (205, 181, 205),
    'thistle4': (139, 123, 139),
    'tomato': (255, 99, 71),
    'tomato1': (255, 99, 71),
    'tomato2': (238, 92, 66),
    'tomato3': (205, 79, 57),
    'tomato4': (139, 54, 38),
    'turquoise': (64, 224, 208),
    'turquoise1': (0, 245, 255),
    'turquoise2': (0, 229, 238),
    'turquoise3': (0, 197, 205),
    'turquoise4': (0, 134, 139),
    'violet': (238, 130, 238),
    'violetred': (208, 32, 144),
    'violetred1': (255, 62, 150),
    'violetred2': (238, 58, 140),
    'violetred3': (205, 50, 120),
    'violetred4': (139, 34, 82),
    'wheat': (245, 222, 179),
    'wheat1': (255, 231, 186),
    'wheat2': (238, 216, 174),
    'wheat3': (205, 186, 150),
    'wheat4': (139, 126, 102),
    'white': (255, 255, 255),
    'whitesmoke': (245, 245, 245),
    'yellow': (255, 255, 0),
    'yellow1': (255, 255, 0),
    'yellow2': (238, 238, 0),
    'yellow3': (205, 205, 0),
    'yellow4': (139, 139, 0),
    'yellowgreen': (154, 205, 50)
}


class TempFileList(object):
    def __init__(self):
        self.flist=[]

    def add(self, fname):
        self.flist.append(fname)

    def cleanup(self):
        for f in self.flist:
            try:
                #print("removing:",f)
                os.remove(f)
            except:
                pass

    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        self.cleanup()
    def __del__(self):
        self.cleanup()

_default_conf={
    'viewer':'eog', # probably available on most systems
}
def _load_config():
    import yaml

    conf=_default_conf
    fname=os.path.expanduser('~/.pyxtools')

    if os.path.exists(fname):
        with open(fname) as fobj:
            tconf=yaml.load(fobj)
        conf.update(tconf)
    return conf

_tflist=TempFileList()
_config=_load_config()
