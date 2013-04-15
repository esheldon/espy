"""
Run the regress.py with type s2n to create the polynomial fits
used to detrend the data.  Then correct using

    g1,g2 = regress.detrend_g_vs_s2n(s2n=, g1=, g2=, 
                                     gmix_run=, camcol=, sratio=)

"""

from . import pipe
from . import files
from . import sweep
from . import collate
from . import regress
from . import cuts
