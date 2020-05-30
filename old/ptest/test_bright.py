import sys
import numpy

# 0 core
# 1 outskirts
# 2 linear structures
# 3 unacceptable

fname = "/Users/esheldon/tmp/test_bright_15_70_15_85.dat"
def stats():

    flags = numpy.fromfile(fname, dtype=numpy.int16, sep='\n')

    print flags

    (wcore,) = numpy.where( (flags & 2**0) != 0 )
    (wskirt,) = numpy.where( (flags & 2**1) != 0 )
    (wlin,) = numpy.where( (flags & 2**2) != 0 )
    (wbad,) = numpy.where( (flags & 2**3) != 0 )
    
    print "Total num: ",flags.size
    print "core: %4.2f%%" % (wcore.size/float(flags.size)*100)
    print "outskirts: %4.2f%%" % (wskirt.size/float(flags.size)*100)
    print "linear: %4.2f%%" % (wlin.size/float(flags.size)*100)
    print "Unacceptable: %4.2f%%" % (wbad.size/float(flags.size)*100)
    
if __name__ == "__main__":
    
    fobj = open(fname, "w")

    print "Enter a comma-separated list of bits"
    print "-----------------------------------------------------"
    input = sys.stdin.readline()

    while input[0] != 'q':

        print "you entered ",input
        bits = numpy.array( eval('['+input+']') )

        if bits[0] == -1:
            print "No Flags: ",0
            flags = numpy.array([0])
        else:
            flags = 2**bits
            print "Flags = ",flags.sum()

        fobj.write("%d\n" % (flags.sum()))

        print "-----------------------------------------------------"
        input = sys.stdin.readline()

    fobj.close()
