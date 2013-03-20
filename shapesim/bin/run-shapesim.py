"""
    %prog run is2 ie/is2n [itrial]

"""

import sys
import os
from esutil.misc import wlog
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)

    itrial=None
    if len(args) > 3:
        itrial = int(args[3])

    run=args[0]
    is2 = int(args[1])
    ie_or_is2n = int(args[2])

    if run[0:5] == 'deswl':
        sim=shapesim.deswl_sim.DESWLSim(run)
    elif run[0:8] == 'gmix-fit':
        sim=shapesim.gmix_fit_sim.GMixFitSim(run)
    elif run[0:4] == 'gmix':
        sim=shapesim.gmix_em_sim.GMixEMSim(run)
    elif 'mixmc' in run or 'bayes' in run or 'mca' in run:
        sim=shapesim.bayesfit_sim.BayesFitSim(run)
        c = shapesim.read_config(run)
        cs = shapesim.read_config(c['sim'])
        if cs['orient'] == 'random':
            if itrial is None:
                raise ValueError("expect itrial for bayes orient=random")
            sim.process_bruteforce_trial_by_s2n(is2, ie_or_is2n, itrial)
            return
    elif 'shd' in run:
        sim=shapesim.shd_sim.SHDirectSim(run)
    elif 'mcm-' in run:
        if itrial is None:
            raise ValueError("send itrial for mcm")
        sim=shapesim.mcm_sim.MCMSim(run)
        sim.process_trial(is2, ie_or_is2n, itrial)

        return

    elif 'stack-' in run:
        if itrial is None:
            raise ValueError("send itrial for stack")

        if 'admom' in run:
            sim=shapesim.stack_sim.AdmomSim(run, is2=is2, is2n=ie_or_is2n, itrial=itrial)
            sim.run()
            sim.write()
        else:
            sim=shapesim.stack_sim.StackSim(run, is2=is2, is2n=ie_or_is2n, itrial=itrial)
            sim.generate_stacks()
            sim.fit_stacks()
            sim.write()
        return

    else:
        raise ValueError("Don't know about run '%s'" % run)

    sim.process_trials(is2, ie_or_is2n, itrial=itrial)


main()
