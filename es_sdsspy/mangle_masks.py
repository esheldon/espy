import os

def load(type, subtype=None, veto=False, verbose=False, code='c'):
    #import mangle
    fname = mask_name(type,subtype=subtype)
    if code == 'wrap':
        import mangle_wrap
        return mangle_wrap.Mangle(fname,veto=veto)
    else:
        import mangle
        return mangle.Mangle(fname,verbose=verbose)

def mask_dir():
    mdir=os.getenv('MASK_DIR')
    if mdir is None:
        raise ValueError("MASK_DIR is not defined")
    mdir = os.path.join(mdir, 'mangle')
    return mdir

def mask_name(type, subtype=None):
    mdir=mask_dir()

    fname=type

    if subtype is not None:
        fname += '-%s' % subtype

    fname += '.ply'

    fname = os.path.join(mdir, fname)
    return fname


