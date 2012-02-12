import sys
import es_sdsspy

if len(sys.argv) < 5:
    print 'usage: generate-randoms.py mapname maptype system nrand > filename'
    print '  system should be eq, gal, or sdss'
    print '  the points to to stdout, so redirect to a file'
    print
    print 'example:'
    print '  generate-randoms.py boss tycho eq 1000000'
    sys.exit(1)

mapname=sys.argv[1]
maptype=sys.argv[2]
system=sys.argv[3]
nrand=int(sys.argv[4])

m=es_sdsspy.stomp_maps.load(mapname,maptype)

carr1,carr2 = m.GenerateRandomPoints(nrand,system)

for c1,c2 in zip(carr1,carr2):
    print '%.16g %.16g' % (c1,c2)
