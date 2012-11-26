import sys
import esutil as eu
from numpy import zeros
import pprint

from optparse import OptionParser
parser=OptionParser(__doc__)

"""
  mem_match_id       >i4  1
  ra                 >f8  239.583329069
  dec                >f8  27.2334129178
  model_mag          >f4  array[5]
  model_magerr       >f4  array[5]
  imag               >f4  13.979
  imag_err           >f4  0.00269815
  zred               >f4  0.0896334
  zred_e             >f4  0.00950781
  bcg_spec_z         >f4  -1.0
  z_spec_init        >f4  0.0805804
  z_init             >f4  0.103065
  z                  >f4  0.0970242
  lambda_chisq       >f4  171.272
  lambda_chisq_e     >f4  3.90487
  r_lambda           >f4  1.11362
  scaleval           >f4  1.01366
  maskfrac           >f4  0.0124726
  c_lambda           >f4  array[4]
  c_lambda_err       >f4  array[4]
  mag_lambda_err     >f4  array[5]
  z_lambda           >f4  0.0946875
  z_lambda_e         >f4  0.00446747
  lnlamlike          >f4  433.908
  lnbcglike          >f4  -16.2419
  lnlike             >f4  417.666
  pzbins             >f4  array[21]
  pz                 >f4  array[21]
  ncross             >i2  1
  chisq              >f4  17.9633
  rmask              >f4  1.67043
  ra_orig            >f8  239.578038802
  dec_orig           >f8  27.2391030452
  w                  >f4  1.05159
  lambda_chisq_c     >f4  171.091
  lambda_chisq_ce    >f4  0.785613
  ncent              >i2  5
  ncent_good         >i2  2
  ra_cent            >f8  array[5]
  dec_cent           >f8  array[5]
  lambda_chisq_cent  
                     >f4  array[5]
  p_bcg              >f4  array[5]
  p_cen              >f4  array[5]
  q_cen              >f4  array[5]
  p_fg               >f4  array[5]
  q_miss             >f4  -1.21325
  p_sat              >f4  array[5]
  p_c                >f4  array[5]
  bcg_ilum           >f4  7.39512
  ilum               >f4  113.608

"""

def make_output(data):
    nokeep=['model_mag','model_magerr','pzbins','pz']

    dt=[]
    for d in data.dtype.descr:
        n=d[0]
        if n in nokeep:
            continue

        if len(data[n].shape) == 1:
            dt += [d]
        else:
            nvals=data[n].shape[1]
            for i in xrange(nvals):
                nout='%s_%d' % (n,i)
                dt += [(nout,d[1])]

    dt += [('umg','f4'),
           ('gmr','f4'),
           ('rmi','f4'),
           ('imz','f4')]


    out=zeros(data.size, dtype=dt)

    for n in data.dtype.names:
        if n in nokeep:
            continue
        if len(data[n].shape) == 1:
            out[n] = data[n]
        else:
            nvals=data[n].shape[1]
            for i in xrange(nvals):
                nout='%s_%d' % (n,i)
                out[nout][:] = data[n][:,i]

    out['umg'] = data['model_mag'][:,0]-data['model_mag'][:,1]
    out['gmr'] = data['model_mag'][:,1]-data['model_mag'][:,2]
    out['rmi'] = data['model_mag'][:,2]-data['model_mag'][:,3]
    out['imz'] = data['model_mag'][:,3]-data['model_mag'][:,4]

    return out

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) == 0:
        sys.stderr.write('tosqlite.py infile outfile\n')
        sys.exit(1)

    infile=args[0]
    dbfile=args[1]

    print 'reading:',infile
    rm=eu.io.read(infile,lower=True)

    eu.numpy_util.to_native(rm,inplace=True)

    out=make_output(rm)

    print 'dbfile:',dbfile
    print 'making table'
    s=eu.sqlite_util.SqliteConnection(dbfile)
    s.array2table(out,'rm')

    s.execute('ALTER TABLE rm ADD comments TEXT')

    print 'adding indices'

    for icol in out.dtype.names:
        print '    ',icol
        s.add_index('rm',icol)

    s.add_index('rm','comments')
           


main()
