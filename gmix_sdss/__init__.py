"""
The order of operations is
    - run on all data
    - sweep results to camcols
    - collate
    - add ra, dec (can occur later)
    - match to photozs (can occur later)
    - add scinv (after matching photoz)
    - do regressions and create detrend files
    - do detrending
    - rotate to equatorial coords

Then one will define source catalog samples in the lensing package
which refer to this catalog.

Details
- Make a config file, e.g. 
        /config/gmix-r02.yaml
- Make the wq scripts 
        /bin/make-wq.py gmix-r02
    - for the bayesian code need to do this by field
    but for LM code can use --bycamcol to have fewer
    jobs

- "sweep" the results
        /bin/make-sweep-wq.py gmix-r02
- collate to columns
        /bin/collate.py gmix-r02
- add ra dec if needed
        /bin/add-radec.py gmix-r02

- Run
        /bin/regress.py
    with type s2n to create the polynomial fits used to
    detrend the data.

- detrend the ellipticities
        /bin/detrend.py gmix-r02

- rotate the detrended ellipticities to the ra,dec frame.  Only the good
    g_dt are used.

        /bin/rotate.py gmix-r02

- match to the photozs (can occur any time after collation)
        bin/match-zphot.py gmix-r02 11

- add scinv values.  Requires photoz match. Takes gmix run and scinv sample

        bin/add-scinv.py gmix-r02 01

    This is super slow


"""

from . import pipe
from . import files
from . import sweep
from . import collate
from . import regress
from . import cuts
