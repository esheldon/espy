import colorsys

def rainbow(num, type='hex'):
    """
    make rainbow colors

    parameters
    ----------
    num: integer
        number of colors
    type: string, optional
        'hex' or 'rgb', default hex
    """
    # not going to 360
    minh = 0.0
    # 270 would go to pure blue
    #maxh = 270.0
    maxh = 285.0

    hstep = (maxh-minh)/(num-1)
    hex_colors=[]
    for i in xrange(num):
        h = minh + i*hstep

        # just change the hue
        r,g,b = colorsys.hsv_to_rgb(h/360.0, 1.0, 1.0)
        r *= 255
        g *= 255
        b *= 255
        if type == 'rgb':
            hex_colors.append(rgb)
        elif type == 'hex':
            hex_colors.append( rgb_to_hex( (r,g,b) ) )
        else:
            raise ValueError("color type should be 'rgb' or 'hex'")

    return hex_colors

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def test_rainbow():
    import numpy
    from biggles import FramedPlot, Points, Curve
    num = 20

    plt = FramedPlot()

    x = numpy.linspace(0.0, 1.0, num)
    y = x**2

    colors = rainbow(num, 'hex')

    for i in xrange(num):
        p = Points([x[i]], [y[i]], type='filled circle', 
                   color=colors[i])
        c = Curve([x[i]],[y[i]], color=colors[i])
        plt.add(p,c)
    plt.show()
