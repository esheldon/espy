#!/usr/bin/env python
import curses
import subprocess

black=0+1
red=1+1
green=2+1
yellow=3+1
blue=4+1
magenta=5+1
cyan=6+1
white=7+1


def qstat(stdscr, username):

    # curses always sets the background to black unless you run this 
    # function.  Default color is -1
    curses.use_default_colors()
    maxy, maxx = stdscr.getmaxyx()
    
    curses.init_pair(1, -1, -1)
    curses.init_pair(blue, curses.COLOR_BLUE, -1)
    curses.init_pair(green, curses.COLOR_GREEN, -1)
    curses.init_pair(yellow, curses.COLOR_YELLOW, -1)
    curses.init_pair(magenta, curses.COLOR_MAGENTA, -1)
    curses.init_pair(red, curses.COLOR_RED, -1)

    x0 = float(maxx)/2
    y0 = float(maxy)/2

    y1 = 0
    x1 = 0
    sizey1 = 3
    win1 = curses.newwin(sizey1, maxx, y1, x1)


    win3 = curses.newwin(maxy-3, maxx, 3, 0)
    maxy3, maxx3 = win3.getmaxyx()

    #center_addstr(win3, 'middle window 3')




    # 30 seconds
    tout_ms = 30*1000
    win1.timeout(tout_ms)
    while 1:
        win3.erase()
        #sout=GetPs()
        sout=GetQstat(username=username)
        sout = sout.split('\n')

        stats = GetQstatStats(sout, username)
        PrintStats(curses, win1, stats, username)
        ii=0
        for line in sout:
            if ii < maxy3 and ii > 1:
                win3.addstr(ii-2, 0, line[0:maxx3])
            ii+=1
        win3.refresh()

        c = win1.getch()
        if c == ord('q'): break  # Exit the while()



def GetQstat(username=None):
    command = 'qstat'
    if username is not None:
        command += ' -u %s' % username
    pobj=subprocess.Popen(command, 
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          shell=True)
    stdout_ret, stderr_ret = pobj.communicate()
    return stdout_ret



def GetQstatStats(sout, username):
    # gather some statistics

    sd={'njobs':0, 
        'nrun': 0, 
        'nq':0,
        'nsusp':0,
        'nexiting':0,
        'nheld':0,
        'nmoving':0,
        'nwait':0}
    for line in sout:
        ls=line.split()

        if len(ls) >= 11:
            if ls[1] == username:
                sd['njobs'] += 1
                jobstat = ls[9]
                if jobstat == "R":
                    sd['nrun'] += 1
                elif jobstat == "Q":
                    sd['nq'] += 1
                elif jobstat == "S":
                    sd['nsusp'] += 1
                elif jobstat == "E":
                    sd['nexiting'] += 1
                elif jobstat == "H":
                    sd['nheld'] += 1
                elif jobstat == "T":
                    sd['nmoving'] += 1
                elif jobstat == "W":
                    sd['nwait'] += 1

    return sd

def PrintStats(curses, win, stats, username):
    win.erase()
    # begin at second line, far left
    win.addstr(1, 0, 'Jobs for %s:' % (username,), curses.color_pair(blue) )
    win.addstr(' %s' % stats['njobs'])

    win.addstr('  Running:', curses.color_pair(green))
    win.addstr(' %s' % stats['nrun'])

    win.addstr('  Queued:', curses.color_pair(yellow))
    win.addstr(' %s' % stats['nq'])

    if stats['nsusp'] > 0:
        win.addstr('  Suspended:', curses.color_pair(magenta))
        win.addstr(' %s' % stats['nsusp'])
    if stats['nexiting'] > 0:
        win.addstr('  Exiting:', curses.color_pair(red))
        win.addstr(' %s' % stats['nexiting'])


    win.refresh()


curses.wrapper(qstat, 'esheldon')
