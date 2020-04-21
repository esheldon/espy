# -*- coding: utf-8 -*-
# Copyright: 2009 Nadia Alramli
# License: BSD
#
# Modification History
#  2012-12-30 Erin Sheldon, Brookhaven National Laboratory
#  - render takes the fraction instead of a percent.
#  - everything is keywords now
#  - Added update() to only render when a change in the percent
#    or message has occurred.
#  2016-11-23
# - made Terminal class from Alramli's terminal module
# - now takes number on construction
# - can use in a with context

"""
Draws an animated terminal progress bar
Usage:
    n = 100
    with ProgressBar(n) as pbar:
        for i in range(n):
            # do something
            p.update()
"""
 
import sys
import time

class Terminal(object):
    COLORS = "BLUE GREEN CYAN RED MAGENTA YELLOW WHITE BLACK".split()

    # List of terminal controls, you can add more to the list.
    CONTROLS = {
        'BOL':'cr', 'UP':'cuu1', 'DOWN':'cud1', 'LEFT':'cub1', 'RIGHT':'cuf1',
        'CLEAR_SCREEN':'clear', 'CLEAR_EOL':'el', 'CLEAR_BOL':'el1',
        'CLEAR_EOS':'ed', 'BOLD':'bold', 'BLINK':'blink', 'DIM':'dim',
        'REVERSE':'rev', 'UNDERLINE':'smul', 'NORMAL':'sgr0',
        'HIDE_CURSOR':'cinvis', 'SHOW_CURSOR':'cnorm'
    }

    # List of numeric capabilities
    VALUES = {
        'COLUMNS':'cols', # Width of the terminal (None for unknown)
        'LINES':'lines',  # Height of the terminal (None for unknown)
        'MAX_COLORS': 'colors',
    }

    def __init__(self):
        self._setup()

    def _setup(self):
        """
        Set the terminal control strings
        """
        import curses

        # Initializing the terminal
        curses.setupterm()

        # Get the color escape sequence template or '' if not supported
        # setab and setaf are for ANSI escape sequences

        bgColorSeq = curses.tigetstr('setab') or curses.tigetstr('setb') or ''
        fgColorSeq = curses.tigetstr('setaf') or curses.tigetstr('setf') or ''

        for color in self.COLORS:
            # Get the color index from curses
            colorIndex = getattr(curses, 'COLOR_%s' % color)

            # Set the color escape sequence after filling the template with
            # index
            setattr(self, color, curses.tparm(fgColorSeq, colorIndex))

            # Set background escape sequence
            setattr(
                self,
                'BG_%s' % color, curses.tparm(bgColorSeq, colorIndex)
            )
        for control in self.CONTROLS:
            # Set the control escape sequence
            setattr(
                self,
                control,
                curses.tigetstr(self.CONTROLS[control]) or '',
            )

        for value in self.VALUES:
            # Set terminal related values
            setattr(self, value, curses.tigetnum(self.VALUES[value]))

class ProgressBar(object):
    """
    Terminal progress bar class
    """
    TEMPLATE = (
    #'%(message)s%(percent)3s%% [%(color)s%(progress)s%(normal)s%(empty)s]\n'
    #PADDING = 8
    '%(message)s%(percent)3s%% ETA: %(eta)8s [%(color)s%(progress)s%(normal)s%(empty)s]\n'
    )
    PADDING = 22
 
    def __init__(self,
                 num,
                 color=None,
                 width=None,
                 block='=',
                 empty=' ',
                 message=None):
        """
        parameters
        ----------
        num: integer
            number of items
        color: string, optional
            terminal color name (BLUE GREEN CYAN RED MAGENTA YELLOW WHITE BLACK)
        width: integer, optional
            width; default is the terminal width.
        block: string, optional
            progress display character (default '=')
        empty: string, optional
            empty bar display character (default ' ')
        """

        self.time_start = time.time()
        self.terminal = Terminal()

        self.num = num

        if color:
            self.color = getattr(self.terminal, color.upper())
        else:
            self.color = ''

        if width and width < self.terminal.COLUMNS - self.PADDING:
            self.width = width
        else:
            # Adjust to the width of the terminal
            self.width = self.terminal.COLUMNS - self.PADDING

        self.block = block
        self.empty = empty
        self.progress = None
        self.lines = 0

        self.current=0
        self.percent_old=-9999
        self.eta_old='blah'
        self.message=message
 
    def update(self):
        """
        print the progress bar if there has been an update.

        An update is when the integer percentage has changed or the message has
        changed.
        """

        self.current += 1

        percent=self.get_percent_message()
        eta=self.get_eta_string()

        if not self.is_updated(percent, eta):
            return

        self.percent_old=percent
        self.eta_old=eta

        self._render(percent, eta)

    def is_updated(self, percent, eta):
        if percent == self.percent_old and eta == self.eta_old:
            return False
        else:
            return True
        return False

    def _render(self, percent, eta):
        """
        Print the progress bar

        parameters
        ----------
        precent: float, optional
            The percent finished.
        """

        inline_msg_len = 0
        if self.message:
            # The length of the first line in the message
            inline_msg_len = len(self.message.splitlines()[0])
            message = '%s ' % self.message
        else:
            message = ''

        if inline_msg_len + self.width + self.PADDING > self.terminal.COLUMNS:
            # The message is too long to fit in one line.
            # Adjust the bar width to fit.
            bar_width = self.terminal.COLUMNS - inline_msg_len -self.PADDING
        else:
            bar_width = self.width
 
        # Check if render is called for the first time
        if self.progress != None:
            self.clear()

        self.progress = (bar_width * percent) / 100

        data = self.TEMPLATE % {
            'percent': percent,
            'eta':eta,
            'color': self.color,
            'progress': self.block * self.progress,
            'normal': self.terminal.NORMAL,
            'empty': self.empty * (bar_width - self.progress),
            'message': message
        }
        sys.stdout.write(data)
        sys.stdout.flush()

        self.data=data
        self.lines = len(data.splitlines())
 
    def get_percent_message(self):
        frac = float(self.current)/self.num
        percent=int(100*frac)
        return percent


    def get_eta_string(self):
        tm = time.time()
        telapse = tm - self.time_start

        time_per = telapse/self.current

        nleft = self.num - self.current

        eta = nleft*time_per

        return time.strftime('%H:%M:%S', time.gmtime(eta))

    def clear(self):
        """Clear all printed lines"""
        sys.stdout.write(
            self.lines * (self.terminal.UP + self.terminal.BOL + self.terminal.CLEAR_EOL)
        )

    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        pass



def _test(message='frac', block='=', empty=' ', color=None):
    n=100

    bar=ProgressBar(
        n,
        block=block,
        empty=empty,
        width=70,
        color=color,
        message=message,
    )
    for i in range(n):
        bar.update()
        time.sleep(0.1)

def test_simple(message='frac'):
    _test(message=message)

def test_utf8(message='fraac'):
    _test(block='▣', empty='□', message=message)

def test_utf8_color(message='frac'):
    _test(block='▣', empty='□', color='green', message=message)


def test_with(message='frac'):

    n=100
    with ProgressBar(n, message='frac') as pbar:
        for i in range(n):
            pbar.update()
            time.sleep(0.1)


