# -*- coding: utf-8 -*-
# Copyright: 2009 Nadia Alramli
# License: BSD
#
# Modification History
#  2012-12-30 Erin Sheldon, Brookhaven National Laboratory
#  - render takes the fraction instead of a percent.
#  - everything is keywords now
#  - Added update() to only render when a change in the percent
#    has occurred.

"""Draws an animated terminal progress bar
Usage:
    p = ProgressBar("blue")
    p.render(percentage, message)
"""
 
import terminal
import sys
 
class ProgressBar(object):
    """Terminal progress bar class"""
    TEMPLATE = (
     '%(percent)-3s%% %(color)s%(progress)s%(normal)s%(empty)s %(message)s\n'
    )
    TEMPLATE_BRACKET = (
     '%(percent)3s%% [%(color)s%(progress)s%(normal)s%(empty)s] %(message)s\n'
    )
    PADDING = 8
 
    def __init__(self, color=None, width=None, block='=', empty=' '):
        """
        color -- color name (BLUE GREEN CYAN RED MAGENTA YELLOW WHITE BLACK)
        width -- bar width (optinal)
        block -- progress display character (default '█')
        empty -- bar display character (default ' ')
        """
        if color:
            self.color = getattr(terminal, color.upper())
        else:
            self.color = ''
        if width and width < terminal.COLUMNS - self.PADDING:
            self.width = width
        else:
            # Adjust to the width of the terminal
            self.width = terminal.COLUMNS - self.PADDING
        self.block = block
        self.empty = empty
        self.progress = None
        self.lines = 0

        self.percent_old=-9999
        self.message_old='nothing'
 
    def update(self, **keys):
        """
        Same as render but only print if there has been an update.

        An update is when the integer percentage has changed or the message has
        changed.

        parameters
        ----------
        frac: float, optional
            The fraction finished.  Percent is 100*frac
        message: string, optional
            A message to print to the right
        """

        percent, message=self.get_percent_message(**keys)

        if not self.is_updated(percent=percent, message=message):
            return
        self.percent_old=percent
        self.message_old=message


        self.render(**keys)

    def is_updated(self, percent=None, message=None):
        if percent is not None:
            if percent == self.percent_old:
                return False
            else:
                return True
        if message is not None:
            if message==self.message_old:
                return False
            else:
                return True
        return False

    def render(self, **keys):
        """
        Print the progress bar

        parameters
        ----------
        frac: float, optional
            The fraction finished.  Percent is 100*frac
        message: string, optional
            A message to print to the right
        """

        percent, message=self.get_percent_message(**keys)

        inline_msg_len = 0
        if message:
            # The length of the first line in the message
            inline_msg_len = len(message.splitlines()[0])
        if inline_msg_len + self.width + self.PADDING > terminal.COLUMNS:
            # The message is too long to fit in one line.
            # Adjust the bar width to fit.
            bar_width = terminal.COLUMNS - inline_msg_len -self.PADDING
        else:
            bar_width = self.width
 
        # Check if render is called for the first time
        if self.progress != None:
            self.clear()
        self.progress = (bar_width * percent) / 100
        data = self.TEMPLATE_BRACKET % {
            'percent': percent,
            'color': self.color,
            'progress': self.block * self.progress,
            'normal': terminal.NORMAL,
            'empty': self.empty * (bar_width - self.progress),
            'message': message
        }
        sys.stdout.write(data)
        sys.stdout.flush()
        # The number of lines printed
        self.lines = len(data.splitlines())
 
    def get_percent_message(self, **keys):
        frac=keys.get('frac',None)
        message=keys.get('message','')

        if frac is not None:
            percent=int(100*frac)
        else:
            percent=None

        return percent, message

    def clear(self):
        """Clear all printed lines"""
        sys.stdout.write(
            self.lines * (terminal.UP + terminal.BOL + terminal.CLEAR_EOL)
        )


