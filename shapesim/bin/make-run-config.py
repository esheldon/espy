"""
    %prog [options] sim_name runtype

"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--s2n-fac',default=None)
parser.add_option('--s2n-vals',default=None)
parser.add_option('--s2n',default=None)
parser.add_option('--ie',default=None,
                  help="ie for bys2n runs")
parser.add_option('--run-name',default=None,
                  help="force run name")
parser.add_option('--overwrite',action='store_true',
                  help="force overwrite")

parser.add_option('--dryrun',action='store_true',
                  help="just print to the screen")




gg_byellip_template="""
run: %(run_name)s
sim: %(sim_name)s

# we will use all the s2 values from the sim, a set of s/n values and a single
# ellip value
runtype: byellip

generic_prior: true

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: false
retrim_fluxfrac: null
s2ncalc_fluxfrac: null

s2n: %(s2n)s
s2n_psf: 1.0e+8

use_cache: false
add_to_cache: false

verbose: false

# we should try with higher values to see what happens
ngauss_psf: 1
ngauss_obj: 1

coellip_psf: true
coellip_obj: true

maxtry: 1
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""

gg_bys2n_template="""
run: %(run_name)s
sim: %(sim_name)s

runtype: bys2n

generic_prior: true

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: false
retrim_fluxfrac: null
s2ncalc_fluxfrac: null

ie: %(ie)s  # %(e_val)0.2f
s2nvals: [%(s2n_vals)s]

s2n_psf: 1.0e+8

use_cache: false
add_to_cache: false

verbose: false

ngauss_psf: 1
ngauss_obj: 1

coellip_psf: true
coellip_obj: true

maxtry: 1
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""




edg_byellip_template="""
run: %(run_name)s
sim: %(sim_name)s

# we will use all the s2 values from the sim, a set of s/n values and a single
# ellip value
runtype: byellip

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: true
retrim_fluxfrac: 0.9973
s2ncalc_fluxfrac: null

s2n: %(s2n)s
s2n_psf: 1.0e+8

use_cache: true
add_to_cache: false

verbose: false

# we should try with higher values to see what happens
ngauss_psf: 2
ngauss_obj: 3

coellip_psf: true
coellip_obj: true

maxtry: 2
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""



edg_bys2n_template="""
run: %(run_name)s
sim: %(sim_name)s

# we will use all the s2 values from the sim, a set of s/n values and a single
# ellip value
runtype: bys2n

generic_prior: true

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: true
retrim_fluxfrac: 0.9973
s2ncalc_fluxfrac: null

ie: %(ie)s  # %(e_val)0.2f
s2nvals: [%(s2n_vals)s]

s2n_psf: 1.0e+8

use_cache: true
add_to_cache: false

verbose: false

# we should try with higher values to see what happens
ngauss_psf: 2
ngauss_obj: 3

coellip_psf: true
coellip_obj: true

maxtry: 2
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""

et_bys2n_template="""
run: %(run_name)s
sim: %(sim_name)s

# we will use all the s2 values from the sim, a set of s/n values and a single
# ellip value
runtype: bys2n

generic_prior: true

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: true
retrim_fluxfrac: 0.9973
s2ncalc_fluxfrac: null

ie: %(ie)s  # %(e_val)0.2f
s2nvals: [%(s2n_vals)s]

s2n_psf: 1.0e+8

use_cache: true
add_to_cache: false

verbose: false

ngauss_psf: 3
ngauss_obj: 3

coellip_psf: true
coellip_obj: true

maxtry: 2
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""

dt_bys2n_template="""
run: %(run_name)s
sim: %(sim_name)s

# we will use all the s2 values from the sim, a set of s/n values and a single
# ellip value
runtype: bys2n

generic_prior: true

s2n_method: admom
s2n_fac: %(s2n_fac)s

retrim: true
retrim_fluxfrac: 0.9973
s2ncalc_fluxfrac: null

ie: %(ie)s  # %(e_val)0.2f
s2nvals: [%(s2n_vals)s]

s2n_psf: 1.0e+8

use_cache: true
add_to_cache: false

verbose: false

ngauss_psf: 3
ngauss_obj: 3

coellip_psf: true
coellip_obj: true

maxtry: 2
maxtry_psf: 1

# this is for the gaussian PSF fit, since admom gets gaussians too perfectly
# trying out doing this automatically as needed
randomize: true

# number of times to retry when a trial fails.  This generates
# a new trial, unlike max_retry above which retries with a 
# randomized guess
itmax: 100

# set to null for new seed each time. Good when adding to the cache
seed: null
"""





shortnames={'turb':'t',
            'gauss':'g',
            'dgauss':'dg',
            'exp':'e',
            'dev':'d'}

def _get_run_name(sim_name, runid):
    simback=sim_name.replace('sim-','')
    run_name = 'gmix-fit-%sr%02d' % (simback,runid)
    return run_name
def get_run_name_and_file(sim_name):

    runid=1
    run_name = _get_run_name(sim_name, runid)
    fname = shapesim.get_config_file(run_name)
    while os.path.exists(fname):
        runid+=1
        run_name = _get_run_name(sim_name, runid)
        fname = shapesim.get_config_file(run_name)

    return run_name, fname

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    sim_name=args[0]
    run_type=args[1]

    cs = shapesim.read_config(sim_name)
    objmodel = cs['objmodel']
    psfmodel = cs['psfmodel']
    simtype = shortnames[objmodel]+shortnames[psfmodel]

    if options.run_name is not None:
        run_name=options.run_name
        fname = shapesim.get_config_file(run_name)
        if os.path.exists(fname) and not options.overwrite:
            raise ValueError("file already exists: %s" % fname)
    else:
        run_name, fname = get_run_name_and_file(sim_name)

    s2n_fac=options.s2n_fac
    if s2n_fac is None:
        raise ValueError("send --s2n-fac")
    if run_type == 'bys2n':
        ie=options.ie
        if ie is None:
            raise ValueError("send --ie for bys2n")

        tmps2,e_val = shapesim.get_s2_e(cs, 0, int(ie))

        #s2n_vals_def='5,10,15,20,25,30,40,50,60,70,80,90,100'
        s2n_vals_def='10,15,20,25,30,40,50,60,70,80,90,100'
        s2n_vals=options.s2n_vals
        if s2n_vals is None:
            s2n_vals=s2n_vals_def
        if simtype == 'edg':
            text=edg_bys2n_template % {'run_name':run_name,
                                       'sim_name':sim_name,
                                       's2n_fac':s2n_fac,
                                       's2n_vals':s2n_vals,
                                       'ie':ie,
                                       'e_val':e_val}
        elif simtype == 'et':
            text=et_bys2n_template % {'run_name':run_name,
                                      'sim_name':sim_name,
                                      's2n_fac':s2n_fac,
                                      's2n_vals':s2n_vals,
                                      'ie':ie,
                                      'e_val':e_val}
        elif simtype == 'dt':
            text=dt_bys2n_template % {'run_name':run_name,
                                      'sim_name':sim_name,
                                      's2n_fac':s2n_fac,
                                      's2n_vals':s2n_vals,
                                      'ie':ie,
                                      'e_val':e_val}
        elif simtype == 'gg':
            text=gg_bys2n_template % {'run_name':run_name,
                                      'sim_name':sim_name,
                                      's2n_fac':s2n_fac,
                                      's2n_vals':s2n_vals,
                                      'ie':ie,
                                      'e_val':e_val}

        else:
            raise ValueError("support other sim types")
    elif run_type == 'byellip':
        s2n=options.s2n
        if s2n is None:
            raise ValueError("send s2n for byellip")
        if simtype == 'edg':
            text=edg_byellip_template
        elif simtype == 'gg':
            text=gg_byellip_template

        else:
            raise ValueError("support other sim types")

        text=text % {'run_name':run_name,
                     'sim_name':sim_name,
                     's2n_fac':s2n_fac,
                     's2n':s2n}

    else:
        raise ValueError("support others?")

    print text
    print 'fname:',fname
    if not options.dryrun:
        with open(fname,'w') as fobj:
            fobj.write(text)
    else:
        print 'dryrun'
main()
