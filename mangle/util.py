import numpy

def construct_cap(num=1):
    return numpy.zeros(num,dtype=[('x','3f8'),('cm','f8')])

def is_cap_used(use_caps, bit):
    return (use_caps & 2**bit) > 0

def set_use_caps(polygon, bitlist, 
                 add=False, 
                 allow_doubles=False, 
                 tol=1.e-10, 
                 allow_neg_doubles=False):

    if not add:
        # these bit flags are not being ored with the existing list, they
        # will just be set
        polygon['use_caps'] = 0

    for bit in bitlist:
        # or with what is already set
        polygon['use_caps'] |= 2**bit

    
    if not allow_doubles:
        for ipoly in polygon.size:
            ncaps = polygon['ncaps'][ipoly]
            for i in range(ncaps):

                if is_cap_used(polygon['use_caps'][ipoly], i):
                    
                    for j in range(i+1, ncaps):

                        if is_cap_used(polygon['use_caps'][ipoly], j):

                            diff2 = (polygon['caps'][ipoly][i]['x'] - polygon['caps'][j]['x'])**2
                            if (diff2.sum() < tol**2):

                                cmdiff = \
                                    polygon['caps'][ipoly][i]['cm'] - polygon['caps'][j]['cm']
                                cmadd = \
                                    polygon['caps'][ipoly][i]['cm'] + polygon['caps'][j]['cm']

                                if (numpy.abs(cmdiff) < tol) \
                                      or ( (numpy.abs(cmadd) < tol) and not allow_neg_doubles ):

                                    # unset the jth bit.  Could also use ~(1 << j)
                                    polygon['use_caps'][ipoly] = \
                                        (polygon['use_caps'][ipoly] & ~(2**j))
