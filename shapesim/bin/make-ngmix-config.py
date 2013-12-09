"""
    %prog sim run fit_model
"""

import sys
import os
import pprint

import shapesim

_template="""
run: "%(run)s"
sim: "%(sim)s"

fit_model: "%(fit_model)s"

nwalkers: 40
burnin:   400
nstep:    200
mca_a:    3.0

# we normalize splits by split for is2n==0
desired_err: %(desired_err)s
nsplit0: %(nsplit0)s

s2n_vals: %(s2n_vals)s 
"""

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--clobber',action='store_true',
                  help="clobber existing file")

def get_pars(fit_model):
    pars={'nsplit0':5000}
    if fit_model=='exp':
        pars['desired_err']='5.0e-05'
        pars['s2n_vals']="[15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400]"
    elif fit_model=='dev':
        pars['desired_err']='1.0e-04'
        pars['s2n_vals']="[15,25,40,60,80,100,140,180,250,500,1000]"
    else:
        raise ValueError("bad fit model: '%s'" % fit_model)

    return pars

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    sim=args[0]
    run=args[1]
    fit_model=args[2]

    out_fname=shapesim.shapesim.get_config_file(run)
    print "writing:",out_fname
    if os.path.exists(out_fname) and not options.clobber:
        raise RuntimeError("file exists: %s" % out_fname)

    pars=get_pars(fit_model)

    pars['sim'] = sim
    pars['run'] = run
    pars['fit_model'] = fit_model

    text=_template % pars

    with open(out_fname,'w') as fobj:
        fobj.write(text)

main()
