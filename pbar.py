"""
modified version of the original simple tqdm from

https://github.com/noamraph/tqdm
"""
__all__ = ['PBar', 'prange']

import sys
import time

def PBar(iterable, desc='', total=None, leave=True, file=sys.stderr,
         mininterval=0.5, miniters=1, n_bars=20):
    """
    Get an iterable object, and return an iterator which acts exactly like the
    iterable, but prints a progress meter and updates it every time a value is
    requested.

    parameters
    ----------
    desc: string, optional
        An optional short string, describing the progress, that is added
        in the beginning of the line.
    total: int, optional
        Optional number of expected iterations. If not given,
        len(iterable) is used if it is defined.
    file: file-like object, optional
        A file-like object to output the progress message to. Default
        stderr
    leave: bool, optional
        If True, leave the remaining text from the progress.  If False,
        delete it.
    mininterval: float, optional
        default 0.5
    miniters: int, optional
        default 1

        If less than mininterval seconds or miniters iterations have passed since
        the last progress meter update, it is not updated again.
    """
    if total is None:
        try:
            total = len(iterable)
        except TypeError:
            total = None
    
    prefix = desc+': ' if desc else ''
    
    sp = StatusPrinter(file)
    sp.print_status(prefix + format_meter(0, total, 0, n_bars=n_bars))
    
    start_t = last_print_t = time.time()
    last_print_n = 0
    n = 0
    for obj in iterable:
        yield obj
        # Now the object was created and processed, so we can print the meter.
        n += 1
        if n - last_print_n >= miniters:
            # We check the counter first, to reduce the overhead of time.time()
            cur_t = time.time()
            if cur_t - last_print_t >= mininterval:
                sp.print_status(prefix + format_meter(n, total, cur_t-start_t, n_bars=n_bars))
                last_print_n = n
                last_print_t = cur_t
    
    if not leave:
        sp.print_status('')
        sys.stdout.write('\r')
    else:
        if last_print_n < n:
            cur_t = time.time()
            sp.print_status(prefix + format_meter(n, total, cur_t-start_t,n_bars=n_bars))
        file.write('\n')


def prange(*args, **kwargs):
    """
    A shortcut for writing PBar(range()) on py3 or Pbar(xrange()) on py2

    e.g.

    import time
    from pbar import prange
    for i in prange(20):
        print(i)
        time.sleep(0.1)
    """
    try:
        f = xrange
    except NameError:
        f = range
    
    return PBar(f(*args), **kwargs)


def format_interval(t):
    mins, s = divmod(int(t), 60)
    h, m = divmod(mins, 60)
    if h:
        return '%d:%02d:%02d' % (h, m, s)
    else:
        return '%02d:%02d' % (m, s)


def format_meter(n, total, elapsed, n_bars=20):
    # n - number of finished iterations
    # total - total number of iterations, or None
    # elapsed - number of seconds passed since start
    if n > total:
        total = None
    
    elapsed_str = format_interval(elapsed)
    #rate = '%5.2f' % (n / elapsed) if elapsed else '?'
    
    if total:
        frac = float(n) / total
        
        bar_length = int(frac*n_bars)
        bar = '#'*bar_length + '-'*(n_bars-bar_length)
        
        percentage = '%3d%%' % (frac * 100)
        
        left_str = format_interval(elapsed / n * (total-n)) if n else '?'
        
        #return '|%s| %d/%d %s [elapsed: %s left: %s, %s iters/sec]' % (
        #    bar, n, total, percentage, elapsed_str, left_str, rate)
        return '|%s| %d/%d %s [elapsed: %s left: %s]' % (
            bar, n, total, percentage, elapsed_str, left_str)
    
    else:
        return '%d [elapsed: %s]' % (n, elapsed_str)


class StatusPrinter(object):
    def __init__(self, file):
        self.file = file
        self.last_printed_len = 0
    
    def print_status(self, s):
        self.file.write('\r'+s+' '*max(self.last_printed_len-len(s), 0))
        self.file.flush()
        self.last_printed_len = len(s)



