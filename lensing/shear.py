from . import util

class Shear:
    """
    A class for shapes and shears

    Examples
    --------
        s1 = Shear()
        s1.set_e1e2(0.2, 0.3)
        
        # read only properties
        print s1.e1, s1.e2, s1.g1, s1.g2

        # This performs proper shear addition
        s3 = s1 + s2
        s4 = s1 - s2
    """
    def __init__(self, e1=None, e2=None, g1=None, g2=None):
        if e1 is not None and e2 is not None:
            self.set_e1e2(e1,e2)
        elif g1 is not None and g2 is not None:
            self.set_g1g2(g1,g2)
        else:
            if e1 is not None or e2 is not None:
                raise ValueError("send both e1,e2")
            elif g1 is not None or g2 is not None:
                raise ValueError("send both g1,g2")
            else:
                self.set_e1e2(0.0,0.0)

    def set_e1e2(self, e1, e2):
        etot = e1**2 + e2**2
        if etot > 1:
            raise ValueError("ellipticity must be <= 1")
        self._e1 = e1
        self._e2 = e2
        self._etot = etot

        g1,g2 = util.e1e2_to_g1g2(self._e1, self._e2)

        self._g1, self._g2 = g1,g2

    def set_g1g2(self, g1, g2):
        e1,e2 = util.g1g2_to_e1e2(g1, g2)
        self.set_e1e2(e1,e2)
        etot = e1**2 + e2**2
        if etot > 1:
            raise ValueError("ellipticity must be <= 1")

        self._g1,self._g2 = g1,g2
        self._e1 = e1
        self._e2 = e2
        self._etot = etot

    def get_e1e2(self):
        return self._e1, self._e2
    def get_e1(self):
        return self._e1
    def get_e2(self):
        return self._e2

    def get_g1g2(self):
        return self._g1, self._g2
    def get_g1(self):
        return self._g1
    def get_g2(self):
        return self._g2

    # read only properties, e.g. shear.e1
    e1 = property(get_e1) 
    e2 = property(get_e2) 
    g1 = property(get_g1) 
    g2 = property(get_g2) 

    def __neg__(self):
        s = Shear()
        s.set_e1e2(-self._e1, -self._e2)
        return s

    def __add__(self, s):
        if not isinstance(s,Shear):
            raise ValueError("Only Shear objects can be added")
        e1o,e2o = util.add_distortion(self._e1, self._e2,
                                      s._e1, s._e2)
        sout = Shear()
        sout.set_e1e2(e1o,e2o)
        return sout
    def __sub__(self, s):
        return self.__add__(-s)
    def __repr__(self):
        return '(%.16g, %.16g)' % (self._e1,self._e2)
