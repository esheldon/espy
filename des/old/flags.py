from esutil.numpy_util import where1

# input columns database, apply "good" selection
def select_good_me_bycol(c):
    flags   = c['shear_flags'][:]
    flagsin = c['input_flags'][:]
    flagswt = c['flags_weight'][:]
    w = where1((flags == 0) & (flagsin == 0) & (flagswt==0))

    return w

def select_good_se_bycol(c):
    flags   = c['shear_flags'][:]
    w = where1(flags == 0)

    return w
