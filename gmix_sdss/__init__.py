"""
- Make a config file, e.g. /config/gmix-r02.yaml
- Make the wq scripts /bin/make-wq.py gmix-r02
    - for the bayesian code need to do this by field
    but for LM code can use --bycamcol to have fewer
    jobs
- "sweep" the results
    /bin/make-sweep-wq.py gmix-r02
- collate to columns
    /bin/collate.py gmix-r02

- Run the /bin/regress.py with type s2n to create the polynomial fits used to
  detrend the data.  Then correct using

    g1,g2 = regress.detrend_g_vs_s2n(s2n=, g1=, g2=, 
                                     gmix_run=, camcol=, sratio=)
- detrend the ellipticities
    /bin/detrend.py

- rotate the detrended ellipticities to the ra,dec frame
    /bin/rotate.py

- match to the photozs
    bin/match-zphot.py gmix-r02 11

(- run /lensing/bin/add-scinv.py from the lensing package, 
    will require a cosmology declared in a config)


- use cuts.py and detrend from regress.py in the scat.py module
    from the lensing package to make a source catalog
"""

from . import pipe
from . import files
from . import sweep
from . import collate
from . import regress
from . import cuts
