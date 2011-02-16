"""
"""

import lensing

class CatalogConverter:
    """
    Read in a catalog and produce an "lcat" binary file for input
    to objshear.  This takes in a class that deals with the reading
    and writing of the files.
    """
    def __init__(self, type):
        if type not in ['scat','lcat']:
            raise ValueError("type must be 'scat' or 'lcat'")
        self.type = type

    def convert(self, objclass, version, sample):
        fname = lensing.files.sample_file(self.type,sample)

        oc = objclass(version)
        oc.create_objshear_input(fname)


