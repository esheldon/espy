# coding: utf-8
# Copyright (c) 2008-2011 Volvox Development Team
#
# Hacked to be a class and tty aware, by Erin Sheldon, BNL
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Author: Konstantin Lepa <konstantin.lepa@gmail.com>

"""ANSII Color formatting for output in terminal."""

from __future__ import print_function
import os,sys


__ALL__ = [ 'Colorizer', 'colored', 'cprint' ]

VERSION = (2, 0, 0)

ATTRIBUTES = dict(
        list(zip([
            'bold',
            'dark',
            '',
            'underline',
            'blink',
            '',
            'reverse',
            'concealed'
            ],
            list(range(1, 9))
            ))
        )
del ATTRIBUTES['']


HIGHLIGHTS = dict(
        list(zip([
            'on_grey',
            'on_red',
            'on_green',
            'on_yellow',
            'on_blue',
            'on_magenta',
            'on_cyan',
            'on_white'
            ],
            list(range(40, 48))
            ))
        )


COLORS = dict(
        list(zip([
            'grey',
            'red',
            'green',
            'yellow',
            'blue',
            'magenta',
            'cyan',
            'white',
            ],
            list(range(30, 38))
            ))
        )


RESET = '\033[0m'

class Colorizer(object):
    def __init__(self, check_tty=False):
        """
        parameters
        ----------
        check_tty: bool
            If True and is not a real terminal (tty) then
            all text is passed through without colors
        """

        self.check_tty=check_tty

        if self.check_tty:
            if sys.stdout.isatty():
                # we are running in a real terminal
                self.docolors=True
            else:
                # we are being piped or redirected
                self.docolors=False
        else:
            self.docolors=True

    def __call__(self, text, color=None, on_color=None, attrs=None):
        """Colorize text.  
        
        If check_tty=True and not a tty, text is passed through unchanged

        Available text colors:
            red, green, yellow, blue, magenta, cyan, white.

        Available text highlights:
            on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.

        Available attributes:
            bold, dark, underline, blink, reverse, concealed.

        Example:
            colored=Colorizer()
            colored('Hello, World!', 'red', 'on_grey', ['blue', 'blink'])
            colored('Hello, World!', 'green')
        """

        if self.docolors:

            fmt_str = '\033[%dm%s'
            if color is not None:
                text = fmt_str % (COLORS[color], text)

            if on_color is not None:
                text = fmt_str % (HIGHLIGHTS[on_color], text)

            if attrs is not None:
                for attr in attrs:
                    text = fmt_str % (ATTRIBUTES[attr], text)

            text += RESET

        return text

    def cprint(self, text, color=None, on_color=None, attrs=None, **kwargs):
        """Print colorize text.

        It accepts arguments of print function.
        """

        ctext = self(text, color=color, on_color=on_color, attrs=attrs)
        print(ctext, **kwargs)


def colored(text, color=None, on_color=None, attrs=None, check_tty=False):
    """Colorize text.  
    
    If check_tty=True and not a tty, text is passed through unchanged

    Available text colors:
        red, green, yellow, blue, magenta, cyan, white.

    Available text highlights:
        on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.

    Available attributes:
        bold, dark, underline, blink, reverse, concealed.

    Example:
        colored=Colorizer()
        colored('Hello, World!', 'red', 'on_grey', ['blue', 'blink'])
        colored('Hello, World!', 'green')
    """

    c=Colorizer(check_tty=check_tty)
    return c(text, color=color, on_color=on_color, attrs=attrs)


def cprint(text, color=None, on_color=None, attrs=None, check_tty=False, **kwargs):
    """Print colorize text.

    It accepts arguments of print function.
    """

    c = Colorizer(check_tty=check_tty)
    c.cprint(text, color=color, on_color=on_color, attrs=attrs, **kwargs)


if __name__ == '__main__':
    print('Current terminal type: %s' % os.getenv('TERM'))
    print('Test basic colors:')

    c=Colorizer()
    c.cprint('Grey color', 'grey', on_color='on_white')
    c.cprint('Red color', 'red')
    c.cprint('Green color', 'green')
    c.cprint('Yellow color', 'yellow')
    c.cprint('Blue color', 'blue')
    c.cprint('Magenta color', 'magenta')
    c.cprint('Cyan color', 'cyan')
    c.cprint('White color', 'white')
    print(('-' * 78))

    print('Test highlights:')
    c.cprint('On grey color', on_color='on_grey')
    c.cprint('On red color', on_color='on_red')
    c.cprint('On green color', color='grey', on_color='on_green')
    c.cprint('On yellow color', color='grey', on_color='on_yellow')
    c.cprint('On blue color', on_color='on_blue')
    c.cprint('On magenta color', on_color='on_magenta')
    c.cprint('On cyan color', color='grey', on_color='on_cyan')
    c.cprint('On white color', color='grey', on_color='on_white')
    print('-' * 78)

    print('Test attributes:')
    c.cprint('Bold grey color', 'grey', attrs=['bold'], on_color='on_white')
    c.cprint('Dark red color', 'red', attrs=['dark'])
    c.cprint('Underline green color', 'green', attrs=['underline'])
    c.cprint('Blink yellow color', 'yellow', attrs=['blink'])
    c.cprint('Reversed blue color', 'blue', attrs=['reverse'])
    c.cprint('Concealed Magenta color', 'magenta', attrs=['concealed'])
    c.cprint('Bold underline reverse cyan color', 'cyan',
            attrs=['bold', 'underline', 'reverse'])
    c.cprint('Dark blink concealed white color', 'white',
            attrs=['dark', 'blink', 'concealed'])
    print(('-' * 78))

    print('Test mixing:')
    c.cprint('Underline red on grey color', 'red', 'on_grey',
            ['underline'])
    c.cprint('Reversed green on grey color', 'green', 'on_grey', ['reverse'])

