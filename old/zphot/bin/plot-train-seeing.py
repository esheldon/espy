"""
    %prog train_sample

"""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-t","--types",default='primus,vvds,zcosmos,deep2',
                  help="types to plot, default %default")
parser.add_option("--yrange",default='0,0.075',
                  help="y range to plot, default %default")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    sample = args[0]
    
    types=options.types.split(',')
    yrange = options.yrange.split(',')
    yrange = [float(v) for v in yrange]

    print 'plotting types:',types
    print 'yrange:        ',yrange

    t = zphot.training.Training(sample)
    t.plot_seeing(types=types,yrange=yrange)

main()
