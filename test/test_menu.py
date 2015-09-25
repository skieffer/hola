#!/usr/bin/python

import sys, os
import subprocess
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

import basic_tools
from hola.hola import HolaConfig
from hola.logging import LogLevel

DIRS = """
orthowontist
sbgn
metro
random
""".split()

def stitchCols(c1, c2, width=54, gutter=5):
    N = max(len(c1), len(c2))
    c1 += ['']*(N-len(c1))
    c2 += ['']*(N-len(c2))
    #fmt = '%%%ds %%s %%s' % width
    #return '\n'.join([
    #    fmt % (l1, ' '*gutter, l2) for l1, l2 in zip(c1, c2)
    #])
    t = ''
    for l1, l2 in zip(c1, c2):
        m = max(0, width - len(l1)) + gutter
        t += l1 + ' '*m + l2 + '\n'
    return t

def buildMenu():
    menutext = ''
    menudict = {}
    col1, col2 = [], []
    col = col1
    n = 0
    for i, d in enumerate(DIRS):
        if i == 3:
            col = col2
        col.extend(['', d])
        graphs = subprocess.check_output('ls graphs/%s/*.gml' % d, shell=True).split('\n')
        for g in graphs:
            if len(g) == 0: continue
            g = g.split('.')[0]
            parts = g.split('/')
            #name = parts[-1].split('.')[0]
            col.append('    (%2d) %s' % (n, parts[-1]))
            filename = '/'.join(parts[-2:])
            menudict[n] = filename
            n += 1
    menutext = stitchCols(col1, col2)
    return (menutext, menudict)

def main():
    menutext, menudict = buildMenu()
    while True:
        print menutext
        c = raw_input('Graph number, or q to quit> ')
        if c == 'q':
            break
        try:
            n = int(c)
            assert n in menudict
        except:
            print 'ERROR: Enter "q" or a graph number'
            continue
        name = menudict[n]
        config = HolaConfig()
        config.useFastSettings(True)
        config.LOG_LEVEL_GENERAL = LogLevel.TIMING
        config.LOG_LEVEL_TREE_PLACEMENT = LogLevel.TIMING
        basic_tools.doHOLA(name, config=config)
        print '%' * 80

if __name__ == '__main__':
    main()

