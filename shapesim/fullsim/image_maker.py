class SimpleCatalogMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create a catalog for the input sim and pointing number

        This is "simple" becuase a constant shear is used, and shapes are drawn
        from a simple shape distribution with random orientations.  All galaxies
        are, for now, the same size.
        """
        pass

class ImageMaker(dict):
    def __init__(self, simname, cat):
        """
        Create a simulated image based on the input catalog

        parameters
        ----------
        simname: string
            The name of the sim.  Corresponds to a config file
            The config controls things like background noise.
        cat: structured array
            The objects to place.  These must have row,col,model,T,e1,e2,counts
        """
        pass
