import sys
import sdsspy

type = sys.argv[1]

ss = sdsspy.pg.SweepStuffer(type)
ss.create_indices()
