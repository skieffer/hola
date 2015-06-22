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

import adaptagrams.adaptagrams as adg
from ortho import Compass
from logging import LogLevel
from constraints import AlignCo

class AlignmentFlags:
    NOTHING = 0
    HALIGN = 1
    VALIGN = 2
    DELIB = 4
    HINFEAS = 8
    VINFEAS = 16

class AlignmentTable:

    def __init__(self, G, ignore=[]):
        """
        :param G: a Graph, for which this table will represent alignment states
        :param ignore: a list of nodes in G which should be ignored

        We maintain a state table in which you can look up a pair of nodes, by their
        IDs, in either order.

        We initialise the table by considering just the NodeConfig object of the
        given graph. (We do not consider the graph's extraPCs list.)
        """
        self.ignoreNodes = ignore
        self.ignoreIDs = [u.ID for u in self.ignoreNodes]
        self.state = {}
        nodes = G.nodes.values()
        # Initialise the table.
        ids = filter(lambda ID: ID not in self.ignoreIDs, [u.ID for u in nodes])
        for id1 in ids:
            d = {}
            for id2 in ids:
                d[id2] = AlignmentFlags.NOTHING
            self.state[id1] = d
        # Now we can add alignments from the graph's NodeConfig with transitive closure.
        for i, u in enumerate(nodes[:-1]):
            if u.ID in self.ignoreIDs: continue
            for v in nodes[i+1:]:
                if v.ID in self.ignoreIDs: continue
                d = G.nodeConf.getDirec(u, v)
                if d in Compass.horizontal:
                    self.addAlignment(u, v, AlignmentFlags.HALIGN)
                elif d in Compass.vertical:
                    self.addAlignment(u, v, AlignmentFlags.VALIGN)

    def areAlignedIDs(self, uID, vID, flag):
        """
        :param uID: a Node ID
        :param vID: a Node ID
        :param flag: HALIGN or VALIGN
        :return: boolean saying if these nodes are aligned in the dimension named
        """
        return self.state[uID][vID] & flag == flag

    def getAlignmentSetIDs(self, u, flag):
        """
        :param u: a Node
        :param flag: HALIGN or VALIGN
        :return: a list of the IDs of all the nodes that are currently
                 aligned with u, including u
        """
        su = self.state[u.ID]
        Au = [u.ID] + [ID for ID in su if su[ID] & flag == flag]
        return Au

    def addAlignment(self, u, v, flag):
        """
        :param u: a Node
        :param v: a Node
        :param flag: HALIGN or VALIGN
        :return: nothing

        We record the alignment between the two nodes, with transitive closure.
        """
        # Get the alignment set IDs:
        Au = self.getAlignmentSetIDs(u, flag)
        Av = self.getAlignmentSetIDs(v, flag)
        # Now record that everything in Au is aligned with everything in Av.
        for id1 in Au:
            for id2 in Av:
                self.state[id1][id2] |= flag
                self.state[id2][id1] |= flag


    def addAlignedNodes(self, nodes, flag):
        """
        :param nodes: a list of Nodes
        :param flag: HALIGN or VALIGN, indicating how the nodes are aligned
        :return: nothing

        We record the alignments, with transitive closure.
        """
        for u, v in zip(nodes[:-1], nodes[1:]):
            self.addAlignment(u, v, flag)

    def noteInfeasibility(self, u, v, flag):
        if flag == AlignmentFlags.HALIGN:
            self.state[u.ID][v.ID] |= AlignmentFlags.HINFEAS
            self.state[v.ID][u.ID] |= AlignmentFlags.HINFEAS
        elif flag == AlignmentFlags.VALIGN:
            self.state[u.ID][v.ID] |= AlignmentFlags.VINFEAS
            self.state[v.ID][u.ID] |= AlignmentFlags.VINFEAS

    def isMarkedInfeasible(self, u, v, flag):
        if flag == AlignmentFlags.HALIGN:
            infeas = AlignmentFlags.HINFEAS
        elif flag == AlignmentFlags.VALIGN:
            infeas = AlignmentFlags.VINFEAS
        return self.state[u.ID][v.ID] & infeas == infeas

class EventKind:
    OPEN = 0
    CLOSE = 1

class Event:

    def __init__(self, node, coord, kind):
        self.node = node
        self.coord = coord
        self.kind = kind

    def __lt__(self, other):
        return self.coord < other.coord

class PotentialAlignment:

    def __init__(self, u, v, flag):
        self.u, self.v = u, v
        self.flag = flag
        self.cost = dist(u, v)
        # For managing linked list:
        self.removed = False
        self.prev = None
        self.next = None

    def __repr__(self):
        return 'PA: %s %d %d' % (
            {AlignmentFlags.HALIGN:'--', AlignmentFlags.VALIGN:'|'}[self.flag],
            self.u.ID, self.v.ID
        )

    def __lt__(self, other):
        return self.cost < other.cost

    def remove(self):
        if not self.removed:
            p, n = self.prev, self.next
            if p is not None:
                p.next = n
            if n is not None:
                n.prev = p
            self.removed = True

    def writeAlignCo(self):
        dim = {AlignmentFlags.HALIGN: adg.YDIM, AlignmentFlags.VALIGN: adg.XDIM}[self.flag]
        ac = AlignCo(dim, shapes=[[self.u, 0], [self.v, 0]])
        return ac

    def getCompassDirec(self):
        if self.flag == AlignmentFlags.HALIGN:
            if self.u.x < self.v.x:
                return Compass.EAST
            else:
                return Compass.WEST
        else:
            if self.u.y < self.v.y:
                return Compass.SOUTH
            else:
                return Compass.NORTH

    def addToNodeConf(self, G):
        """
        :param G: a Graph
        :return: nothing

        We add this alignment to the graph's node config.
        """
        d = self.getCompassDirec()
        G.nodeConf.setDirec(self.u, self.v, d)

    def addToTable(self, atab):
        """
        :param atab: an AlignmentTable
        :return: nothing

        We add this alignment to the table.
        """
        atab.addAlignment(self.u, self.v, self.flag)

    def noteInfeasibility(self, atab):
        atab.noteInfeasibility(self.u, self.v, self.flag)


def dist(u, v):
    "Manhattan metric"
    ux, uy = u.centre()
    vx, vy = v.centre()
    return abs(ux - vx) + abs(uy - vy)


def doNearAlignments(avgdim, logger, G, atab, config, ignore=[], reattempt=False):
    """
    :param G: a Graph
    :param atab: an AlignmentTable for the graph
    :param config: HolaConfig object
    :param ignore: a list of Nodes to be ignored
    :return: nothing

    We look for nodes that are nearly aligned, and try to align them.
    """
    XDIM, YDIM = adg.XDIM, adg.YDIM
    active = {XDIM: True, YDIM: True}
    nextDim = {XDIM: YDIM, YDIM: XDIM}
    d2a = {XDIM: AlignmentFlags.VALIGN, YDIM: AlignmentFlags.HALIGN}
    dim = XDIM
    nodes = {}
    for node in G.nodes.values():
        if not node in ignore:
            nodes[node.ID] = node
    hkw = config.getKinkWidth(avgdim) / 2.0
    asc = config.getAlignmentScope(avgdim)
    def inScope(u, v):
        return abs(scopeCoord(u) - scopeCoord(v)) <= asc
    while active[XDIM] or active[YDIM]:
        dim = nextDim[dim]
        aflag = d2a[dim]
        active[dim] = False
        # Map from node ID to list of nodes that are candidates for
        # alignment in the current dimension:
        candidates = {}
        # Prepare a list of open and close events for the nodes.
        evts = []
        for node in nodes.values():
            x, X, y, Y = node.boundingBoxxXyY()
            cx, cy = (x+X)/2.0, (y+Y)/2.0
            Ix, Iy = (cx - hkw, cx + hkw), (cy - hkw, cy + hkw)
            z, Z = Ix if dim == XDIM else Iy
            # Old method:
            #z, Z = (x, X) if dim == XDIM else (y, Y)
            evts.append(Event(node, z, EventKind.OPEN))
            evts.append(Event(node, Z, EventKind.CLOSE))
        evts.sort()
        scopeCoord = (lambda u: u.y) if dim == XDIM else (lambda u: u.x)
        # Maintain a list of "open" node IDs:
        openIDs = []
        # Scan through the events:
        for evt in evts:
            if evt.kind == EventKind.OPEN:
                u = evt.node
                uc = candidates.get(u.ID, [])
                for vID in openIDs:
                    v = nodes[vID]
                    if ((not atab.areAlignedIDs(u.ID, vID, aflag))
                        and inScope(u, v)
                    ):
                        vc = candidates.get(vID, [])
                        # We don't want to record the candidate twice, so
                        # just store it in the list for the node with smaller ID.
                        if u.ID < vID:
                            uc.append(v)
                        else:
                            vc.append(u)
                        candidates[vID] = vc
                candidates[u.ID] = uc
                openIDs.append(u.ID)
            else:
                assert(evt.kind == EventKind.CLOSE)
                openIDs.remove(evt.node.ID)
        # Now build a list of PotentialAlignments, one for each candidate.
        pas = []
        for uID, C in candidates.items():
            u = nodes[uID]
            for v in C:
                if reattempt or not atab.isMarkedInfeasible(u, v, aflag):
                    pas.append(PotentialAlignment(u, v, aflag))
        # If there aren't any, then move on.
        if len(pas) == 0: continue
        # Sort by cost, lowest first.
        pas.sort()
        # Build a linked list.
        first = pas[0]
        ext = [None] + pas + [None]
        for i in range(len(pas)):
            o, p, q = ext[i:i+3]
            p.prev, p.next = o, q
        # Build an index, in which a node ID points to a dictionary whose values
        # are all the PAs in which that node is involved. The keys are the IDs of
        # the /other/ node in the PA in question.
        index = {}
        for pa in pas:
            for t, w in [(pa.u, pa.v), (pa.v, pa.u)]:
                d = index.get(t.ID, {})
                d[w.ID] = pa
                index[t.ID] = d
        # Now we choose the "prime" alignments, so called because the process by which
        # we choose them is much like the Sieve of Eratosthenes.
        primePAs = []
        while first is not None:
            # Accept the first PA.
            primePAs.append(first)
            # Let u, v be the two nodes involved.
            u, v = first.u, first.v
            # Cross off any other PA in which u or v is involved.
            for w in [u, v]:
                for pa in index[w.ID].values():
                    pa.remove()
            # Furthermore, for any node w already aligned with u, remove any potential
            # alignment of w with v. Likewise for nodes already aligned with v.
            for t, w in [(u, v), (v, u)]:
                At = atab.getAlignmentSetIDs(t, aflag)
                for ID in At:
                    d = index.get(ID, None)
                    if d is None: continue
                    pa = d.get(v.ID, None)
                    if pa is not None:
                        pa.remove()
            # Advance the pointer to the first PA.
            first = first.next
        # Attempt to apply alignments.
        for pa in primePAs:
            if logger.level >= LogLevel.PROGRESS:
                print pa
            ac = pa.writeAlignCo()
            result = G.project(logger, dim, op=True, pcs=[ac], solidEdges=True)
            if logger.level >= LogLevel.PROGRESS:
                print '    %d' % result
            if result == 0:
                active[dim] = True
                pa.addToNodeConf(G)
                pa.addToTable(atab)
            else:
                pa.noteInfeasibility(atab)
