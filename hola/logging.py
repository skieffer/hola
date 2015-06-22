######################################################################
#
# HOLA == Human-like Orthogonal Layout Algorithm
# This file is part of HOLA.
#
# Copyright (C) 2014-2015  Steve Kieffer
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.
#
#
# Author: Steve Kieffer  <http://skieffer.info>
#
######################################################################

import os
import time

class LogLevel:
    NONE = 0x0
    TIMING = 0x10
    PROGRESS = 0x18
    NODE_IDS_AS_LABELS = 0x20
    STAGE_GRAPHS = 0x30
    FINER_STAGE_GRAPHS = 0x40

    DEBUG = 0x1000


class Logger:

    def __init__(self, folderName, basicFileName):
        self.folderName = folderName
        self.basicFileName = basicFileName
        if len(folderName) > 0:
            self.createFolder(self.folderName)
        self.lastSuffix = ''
        self.graph = None
        self.rootTimer = Timer('ROOT', None)
        self.timerStack = []
        self.numProjections = 0
        self.level = 0
        self.levelStack = []

    def pushLevel(self):
        self.levelStack.append(self.level)

    def popLevel(self):
        self.level = self.levelStack.pop()

    def addProj(self):
        self.numProjections += 1

    def startNewTimer(self, name):
        parent = self.rootTimer if len(self.timerStack) == 0 else self.timerStack[-1]
        timer = Timer(name, parent)
        self.timerStack.append(timer)
        timer.start()

    def stopLastTimer(self):
        self.timerStack[-1].stop()
        self.timerStack.pop()

    # Convenience method:
    def startT(self, name):
        if self.level >= LogLevel.TIMING:
            self.startNewTimer(name)

    # Convenience method:
    def stopT(self):
        if self.level >= LogLevel.TIMING:
            self.stopLastTimer()

    def createFolder(self, fn):
        path = fn
        if not os.path.exists(path):
            os.makedirs(path)

    def writeGML(self, suffix, graph=None, highlightNodes=None):
        if graph is None:
            graph = self.graph
        if graph is None:
            print 'Debug output manager has no graph!'
            return
        path = '%s/%s%s.gml' % (
            self.folderName, self.basicFileName, suffix
        )
        f = open(path,'w')
        f.write(graph.writeGML(highlightNodes=highlightNodes))
        f.close()
        self.lastSuffix = suffix

    def writeSVG(self, suffix, w, h, graph=None):
        if graph is None:
            graph = self.graph
        if graph is None:
            print 'Debug output manager has no graph!'
            return
        path = '%s/%s%s.svg' % (
            self.folderName, self.basicFileName, suffix
        )
        f = open(path,'w')
        f.write(graph.writeSVG(w, h))
        f.close()
        self.lastSuffix = suffix

    def writeRouterSVG(self, router, suffix):
        path = '%s/%s%s_route' % (
            self.folderName, self.basicFileName, suffix
        )
        router.outputInstanceToSVG(path)


class Timer:

    def __init__(self, name, parent):
        self.name = name
        self.parent = parent
        self.startTime = 0
        self.stopTime = 0
        self.kids = []
        if self.parent is not None:
            self.parent.addTimer(self)

    def results(self, pad=''):
        s = ''
        s += pad + 'Timer: %s\n' % self.name
        dur = self.duration()
        s += pad + '  ' + 'Total time: %.5f\n' % dur
        if len(self.kids) > 0:
            s += pad + '  ' + 'Subtotals:\n'
            sts = self.subtotals()
            a = 0
            for st in sts:
                v = sts[st]
                a += v
                s += pad + '    %s: %.5f\n' % (st, v)
            if self.name != 'ROOT':
                s += pad + '    other: %.5f\n' % (dur - a)
            s += pad + '  ' + 'Children:\n'
            for kid in self.kids:
                s += kid.results(pad=pad+'    ')
        return s

    def addTimer(self, kid):
        self.kids.append(kid)

    def start(self):
        self.startTime = time.time()

    def stop(self):
        self.stopTime = time.time()

    def ancestors(self):
        a = self.name
        if self.parent is not None:
            a += ', ' + self.parent.ancestors()
        return a

    def duration(self):
        dur = self.stopTime - self.startTime
        if dur < 0:
            print 'Warning: negative duration %.5f for %s' % (dur, self.ancestors())
        return dur

    def subtotals(self):
        sts = {}
        for kid in self.kids:
            st = sts.get(kid.name, 0)
            st += kid.duration()
            sts[kid.name] = st
        return sts

    def subtotalsRec(self, sts={}):
        for kid in self.kids:
            st = sts.get(kid.name, 0)
            st += kid.duration()
            sts[kid.name] = st
        return sts

    def getTimingDataRec(self):
        dur = self.duration()
        td = TimingData(dur)
        cats = {}
        for kid in self.kids:
            td0 = cats.get(kid.name, TimingData(0))
            td1 = kid.getTimingDataRec()
            cats[kid.name] = td0 + td1
        if len(cats.keys()) > 0 and dur > 0:
            other = dur - sum([d.total for d in cats.values()])
            cats['other'] = TimingData(other)
        td.categories = cats
        return td


class TimingData:

    def __init__(self, total):
        self.total = total
        # Dictionary of the form key:value = string:TimingData
        self.categories = {}

    def write(self, div=1, pad='', grandTotal=0):
        div = float(div)
        s = ''
        if self.total > 0:
            s += '%.5f\n' % (self.total/div)
        cats = self.categories.keys()
        if len(cats) > 0:
            # Move 'other' to the end of the list, if it's in there.
            try:
                cats.remove('other')
            except:
                pass
            else:
                cats.append('other')
            # Write the data.
            for cat in cats:
                td = self.categories[cat]
                s += pad + td.percentageString(self.total) + cat + ': '
                s += td.write(div=div, pad=pad+'    ', grandTotal=self.total)
        return s

    def percentageString(self, grandTotal):
        return '' if grandTotal==0 else ('(%2d%%) ' % int(round(100*self.total/grandTotal)))

    def setCategory(self, name, td):
        self.categories[name] = td

    def __add__(self, other):
        t = self.total + other.total
        S = set(self.categories.keys())
        O = set(other.categories.keys())
        A = S.union(O)
        td = TimingData(t)
        for cat in A:
            s = self.categories.get(cat, TimingData(0))
            o = other.categories.get(cat, TimingData(0))
            td.setCategory(cat, s + o)
        return td

