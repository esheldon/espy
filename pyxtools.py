from pyx import *

def test_basic():

    # can define by hex string or rgb value, normalized to 1
    #yellow=color.rgbfromhexstring('#FFFF00')
    yellow=color.rgb(r=1.0, g=1.0, b=0.0)

    # need to specify one or both of width, height
    key=graph.key.key(pos="br", dist=0.1)
    g = graph.graphxy(
        width=8,
        height=8,
        key=key
    )

    # either provide lists of the individual coordinates
    x=list(range(10))
    y=[xx**2 for xx in x]
    values=graph.data.values(x=x, y=y,title='points')
    symbol1=graph.style.symbol(
        symbol=graph.style._circlesymbol,
        symbolattrs=[deco.filled([yellow])],
    )
    symbol2=graph.style.symbol(
        symbol=graph.style._circlesymbol,
        symbolattrs=[color.rgb.red],
    )

    g.plot(values,[symbol1,symbol2])

    # or provide one list containing the whole points
    #values2=graph.data.points(list(zip(range(10), range(10))), x=1, y=2)
    #g.plot(values2)

    # also a function
    c=graph.data.function("y(x)=x**2", min=0, max=9,title=r'$y=x^2$')
    g.plot(c)

    g.writeEPSfile("points-and-line")
    g.writePDFfile("points-and-line")
    #g.writeSVGfile("points")
