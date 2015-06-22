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

from ortho import Compass
from collections import deque as Deque
import math

class LinkShape:
    r = 0  # looks like Latin lowercase 'r'
    u = 1  # looks like Hangul 'uh' character
    n = 2  # looks like Hangul n
    g = 3  # looks like Hangul g
    i = 4  # looks like Hangul i ("ee" sound)
    j = 5  # looks like Latin 'J' (sans serif)

    bent = [0, 2, 3, 5]
    bentCW = [0, 3, 5, 2]
    straight = [1, 4]

    @classmethod
    def cwBendsFrom(cls, firstBend):
        k0 = cls.bentCW.index(firstBend)
        return cls.bentCW[k0:] + cls.bentCW[:k0]

def applyBendToDir(b, d):
    """
    :param b: a bent LinkShape
    :param d: a cardinal Compass direction
    :return: the cardinal direction you would be going if you came into bend b
             going in direction d, and then followed the bend; None if b and d
             are incompatible
    """
    E, S, W, N = Compass.cwCards
    r, g, j, n = LinkShape.bentCW
    lookup = {
        r: {N: E, W: S},
        g: {E: S, N: W},
        j: {S: W, E: N},
        n: {W: N, S: E}
    }
    return lookup[b].get(d, None)

def cwIncomingDirForBend(b):
    """
    :param b: a bent LinkShape
    :return: the clockwise incoming Compass direction for the given bend type
    """
    E, S, W, N = Compass.cwCards
    r, g, j, n = LinkShape.bentCW
    return {
        r: N, g: E, j: S, n: W
    }[b]


# We have a lookup table for the bend sequences that can occur on bend-minimal
# paths from node A to node Z, leaving A in direction d0 and entering Z in
# direction d1.
#
# We only write the table for the case where direction d1 is E.
#
# To use the table, you must first rotate your situtation so that d1 is E.
# Then compute the relative position of node A w.r.t. node Z. This relative
# position is encoded as an integer p in 0, 1, 2, ..., 8, where these integers
# represent the positions in the following grid:
#            ___ ___ ___
#           | 0 | 1 | 2 |
#           `--- --- ---'
#           | 3 | 4 | 5 |
#           `--- --- ---'
#           | 6 | 7 | 8 |
#           `--- --- ---'
# The idea is that Z sits at position 4, and p will indicate at which position A
# sits, where
#
#     p == 0 means A.x < Z.x ^ A.y < Z.y,
#     p == 1 means A.x == Z.x ^ A.y < Z.y,
#     ...
#     p == 8 means A.x > Z.x ^ A.y > Z.y.
#
# HOWEVER: It turns out the the bend sequences for A at position 1 are identical
# with those for A at position 2; likewise pos 4 is same as pos 5, and pos 7 is
# the same as pos 8.
# Therefore the lookup table in fact has no entries for positions 1, 4, 7.
# For these you must look up 2, 5, 8, respectively.
#
# So with p chosen from {0, 3, 6, 2, 5, 8}, you grab the list of possible bend
# sequences:
#
#     eastFinalLookup[p][d0]
#
# And then you must "rotate" the entries in these sequences back to your original
# orientation, if your d1 was not actually E.
#
# For convenience we provide a function
#
#     lookupMinimalBendSeqs(A, d0, Z, d1)
#
# that does all of the rotation and lookup work for you. You just pass
# it A, d0, Z, d1, and it returns the list of bend sequences for your actual
# situation. See below.
#
E, S, W, N = Compass.cwCards
r, g, j, n = LinkShape.bentCW
eastFinalLookup = {
    0: {
        E: [[g,n]],
        S: [[n]],
        W: [[r,n]],
        N: [[r,g,n], [g,r,n]]
    },
    2: {
        E: [[j,g,r,n], [g,j,r,n], [g,j,n,r]],
        S: [[j,r,n], [j,n,r]],
        W: [[r,n]],
        N: [[g,r,n]]
    },
    3: {
        E: [[]],
        S: [[n,j,r]],
        W: [[n,r,g,n], [r,n,j,r]],
        N: [[r,g,n]]
    },
    5: {
        E: [[j,g,r,n], [g,j,n,r]],
        S: [[j,n,r]],
        W: [[n,g,r,n], [r,j,n,r]],
        N: [[g,r,n]]
    },
    6: {
        E: [[j,r]],
        S: [[n,j,r], [j,n,r]],
        W: [[n,r]],
        N: [[r]]
    },
    8: {
        E: [[g,j,n,r], [j,g,n,r], [j,g,r,n]],
        S: [[j,n,r]],
        W: [[n,r]],
        N: [[g,n,r], [g,r,n]]
    }
}
def lookupMinimalBendSeqs(A, d0, Z, d1):
    """
    :param A: Node at beginning of path
    :param d0: Compass direction in which to depart from A
    :param Z: Node at end of path
    :param d1: Compass direction in which to enter Z
    :return: a list [s0, s1, ..., sk-1] where each si is a list of LinkShapes.
             These are all and only the sequences of bends that can occur on
             a bend-minimal orthogonal route from node A to node Z, with the
             prescribed departure and arrival directions d0 and d1.

             You always get at least one si, but it itself may be empty (meaning
             that the best path has zero bends).
    """
    E, S, W, N = Compass.cwCards
    r, g, j, n = LinkShape.bentCW
    pMap = {
        E: [0, 2, 2, 3, 5, 5, 6, 8, 8],
        S: [6, 3, 0, 8, 5, 2, 8, 5, 2],
        W: [8, 8, 6, 5, 5, 3, 2, 2, 0],
        N: [2, 5, 8, 2, 5, 8, 0, 3, 6]
    }[d1]
    d0Map = {
        E: [0, 1, 2, 3],
        S: [3, 0, 1, 2],
        W: [2, 3, 0, 1],
        N: [1, 2, 3, 0]
    }[d1]
    bendMap = {
        E: {r:r, g:g, j:j, n:n},
        S: {r:g, g:j, j:n, n:r},
        W: {r:j, g:n, j:r, n:g},
        N: {r:n, g:r, j:g, n:j}
    }[d1]
    y = 0 if A.y < Z.y else (1 if A.y == Z.y else 2)
    x = 0 if A.x < Z.x else (1 if A.x == Z.x else 2)
    p = 3*y + x
    return [
        [ bendMap[b] for b in seq ]
        for seq in eastFinalLookup[pMap[p]][d0Map[d0]]
    ]

def shapeOfLink(link):
    """
    :param link: a Node of degree 2
    :return: the LinkShape for the shape of this link
    """
    d = []
    for edge in link.edges.values():
        v = edge.otherEnd(link)
        d.append( Compass.cardinalDirection(link, v) )
    d.sort()
    return {
        '01': LinkShape.r, '02': LinkShape.u, '03': LinkShape.n,
        '12': LinkShape.g, '13': LinkShape.i, '23': LinkShape.j
    }.get('%d%d'%tuple(d), None)

class BendSequence:
    """
    A data structure for managing sequences of bend types, points at which these
    bends should occur (in a given Chain), cost of such a sequence of bends
    (for a given Chain), and incoming and outgoing Compass directions, for non-cycles
    """

    def __init__(self, bendtypes, dIn=None, dOut=None):
        self.bendtypes = bendtypes
        self.bendpoints = []
        self.cost = 0
        self.incomingDirec = dIn
        self.outgoingDirec = dOut

    def __repr__(self):
        s = str(self.bendtypes)
        if self.incomingDirec is not None:
            s += ' entering %d' % self.incomingDirec
        if self.outgoingDirec is not None:
            s += ' exiting %d' % self.outgoingDirec
        return s

class AestheticBend:

    def __init__(self, bendNode, nbrNode1, nbrNode2):
        self.bendNode = bendNode
        self.nbrNode1 = nbrNode1
        self.nbrNode2 = nbrNode2

    def addRoutePointToEdgeInGraph(self, G):
        """
        :param G: a Graph object
        :return: nothing

        Get the Edge in the given graph between the two neighbour nodes of
        this bend. Add to it a route point at the position of the bend point.
        """
        edge = G.getEdgeBtwNodes(self.nbrNode1, self.nbrNode2)
        assert(edge is not None)
        edge.addBendPt(self.bendNode.centre())


class Chain:

    def __init__(self, G, nodes, cycle=False):
        """
        :param G: the Graph to which the Chain belongs
        :param nodes: a list(-like object) of the nodes belonging to the chain, in order
        :param cycle: boolean saying if this chain forms a cycle
        """
        self.graph = G
        self.aestheticBends = []
        if len(nodes) == 0:
            raise Exception("Cannot have an empty chain.")
        self.nodes = list(nodes)
        self.cycle = cycle
        if self.cycle:
            assert(len(nodes) >= 3)
            # For cycles, we always store the nodes in clockwise order.
            # Start by getting the index of a node of minimal y-coord.
            i1 = sorted(list(enumerate(self.nodes)),key=lambda p:p[1].y)[0][0]
            i0 = (i1 - 1) % len(self.nodes)
            i2 = (i1 + 1) % len(self.nodes)
            n0, n1, n2 = map(lambda i: self.nodes[i], [i0, i1, i2])
            if n0.x < n1.x:
                # already clockwise
                pass
            elif n0.x > n1.x:
                # anticlockwise
                self.nodes.reverse()
            else:
                # Part and parcel of the assumption that the cycle even /has/ an interior
                # is the assumption that it is not self-intersecting. Therefore, since
                # both neighbouring nodes n0 and n2 have y-coord >= that of n1, they cannot
                # both have the same x-coord as n1. Therefore...
                assert(n2.x != n1.x)
                if n2.x > n1.x:
                    # already clockwise
                    pass
                else:
                    # anticlockwise
                    self.nodes.reverse()
        # Compute and store the shape of each link.
        self.shapes = map(shapeOfLink, nodes)
        # Determine the sequence of internal edges, as well as the
        # anchor nodes and edges if it is not a cycle,
        # or the "return edge" if it is a cycle.
        self.edges = []
        self.anchorNodeLeft = None
        self.anchorEdgeLeft = None
        self.anchorNodeRight = None
        self.anchorEdgeRight = None
        self.returnEdge = None
        n0 = self.nodes[0]
        e1, e2 = n0.edges.values()
        n1, n2 = e1.otherEnd(n0), e2.otherEnd(n0)
        if len(self.nodes) == 1:
            assert(not self.cycle)
            # In this case 'left' and 'right' are meaningless, so record
            # the nodes and edges in any way.
            self.anchorNodeLeft = n1
            self.anchorEdgeLeft = e1
            self.anchorNodeRight = n2
            self.anchorEdgeRight = e2
        else:
            if n1 == self.nodes[1]:
                self.anchorNodeLeft = n2
                self.anchorEdgeLeft = e2
                e0 = e1
            else:
                self.anchorNodeLeft = n1
                self.anchorEdgeLeft = e1
                e0 = e2
            for n0 in self.nodes[1:]:
                # Append the edge e0 that leads into n0 from the left.
                self.edges.append(e0)
                # And get the next edge.
                E = set(n0.edges.values())
                E.remove(e0)
                for e0 in E: break
            self.anchorEdgeRight = e0
            self.anchorNodeRight = e0.otherEnd(n0)
            if self.cycle:
                # In the case of a cycle, the "anchors" are meaningless, but harmless.
                self.returnEdge = self.anchorEdgeRight

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        s = 'Chain: %s' % [node.ID for node in self.nodes]
        #if self.cycle:
        #    s += ' (cycle)'
        #elif self.isEll():
        #    s += ' (ell)'
        verbose = True
        if verbose:
            s += '\n    Internal edges:\n'
            for edge in self.edges:
                s += '    %s\n' % edge
            if self.cycle:
                s += '    Cycle return edge: %s\n' % self.returnEdge
            else:
                s += '    Left anchor: %s, %s\n' % (
                    self.anchorNodeLeft.ID, self.anchorEdgeLeft
                )
                s += '    Right anchor: %s, %s\n' % (
                    self.anchorNodeRight.ID, self.anchorEdgeRight
                )
        return s

    def addRoutePointsInGraph(self, G):
        for a in self.aestheticBends:
            a.addRoutePointToEdgeInGraph(G)

    def getNode(self, i):
        """
        Together with the getEdge function, this function allows us to have the indices
            0, 1, 2, 3, ...
        refer to the first node in the chain, then the first edge, next node, next edge, ...

        :param i: an even integer from -2 to 2n, where n is the number of nodes in
                  this chain.
        :return: left anchor node for i == -2, self.nodes[i/2] for i from 0 to 2n-2,
                 and right anchor node for i == 2n
        """
        n = len(self.nodes)
        assert(i%2==0 and -2 <= i and i <= 2*n)
        if i == -2:
            return self.anchorNodeLeft
        elif i == 2*n:
            return self.anchorNodeRight
        else:
            return self.nodes[i/2]

    def getEdge(self, i):
        """
        Together with the getNode function, this function allows us to have the indices
            0, 1, 2, 3, ...
        refer to the first node in the chain, then the first edge, next node, next edge, ...

        :param i: an odd integer from -1 to 2n-1, where n is the number of nodes in
                  this chain
        :return: left anchor edge for i == -1, self.edges[(i-1)/2] for i from 1 to 2n-3,
                 and right anchor edge for i == 2n-1
        """
        n = len(self.nodes)
        assert(i%2==1 and -1 <= i and i <= 2*n-1)
        if i == -1:
            return self.anchorEdgeLeft
        elif i == 2*n - 1:
            return self.anchorEdgeRight
        else:
            return self.edges[(i-1)/2]


    def numBends(self):
        return len(filter(lambda sh: sh in LinkShape.bent, self.shapes))

    def isEll(self):
        """
        "Ell-chains" are those with precisely one bent link.
        :return: boolean
        """
        return self.numBends() == 1

    def getNodePairsForRange(self, i0, i1):
        """
        :param i0: index in range [0, len(self.nodes) - 2]
        :param i1: index in range [1, len(self.nodes) - 1]
        :return: node pairs (i0, i0+1), (i0+1, i0+2), ..., (i1-1, i1)
        """
        return [
            (self.nodes[a], self.nodes[a+1])
            for a in range(i0, i1)
        ]

    def nextBendIndex(self, i0):
        """
        :param i0: an index
        :return: the index of the next bent link on or after index i0, or -1 if none
        """
        i1 = -1
        for i, shape in enumerate(self.shapes[i0:]):
            if shape in LinkShape.bent:
                i1 = i
                break
        return i1

    def firstBendIndex(self):
        """
        :return: index of first bent link, or -1 if none
        """
        return self.nextBendIndex(0)

    def bendCost(self, bendtype, i0):
        """
        :param bendtype: a bent LinkShape
        :param i0: a position in the chain -- evens for nodes, odds for edges
        :return: the cost of creating that bend shape at that position, given
                 current geometry.
                 If this Chain is a cycle, then the cost takes into account that
                 the nodes are in clockwise order.
        """
        # First compute the angle alpha for position i0.
        # This is the atan2 for a vector z from point p to point q, where
        # if i0 is an edge then p is centre of node i0 - 1 and q centre of node i0 + 1;
        # if i0 is a node then p and q are points on edges i0 - 1 and i0 + 1 resp, a
        # each a unit distance from centre of node i0.
        if i0 % 2 == 1:
            u, w = self.getNode(i0 - 1), self.getNode(i0 + 1)
            p, q = (u.x, u.y), (w.x, w.y)
        else:
            u, v, w = [self.getNode(i0 + di) for di in [-2, 0, 2]]
            p, q = (u.x - v.x, u.y - v.y), (w.x - v.x, w.y - v.y)
            lp, lq = math.sqrt(p[0]**2 + p[1]**2), math.sqrt(q[0]**2 + q[1]**2)
            p = [c/lp for c in p]
            q = [c/lq for c in q]
        z = [q[0] - p[0], q[1] - p[1]]
        # Get angle in degrees.
        alpha0 = math.atan2(z[1], z[0]) * (180 / math.pi)
        r, g, j, n = LinkShape.bentCW
        if self.cycle:
            # For a cycle each type of bend has a specific angle associated with it,
            # so you can be up to +/-180 degrees off.
            beta = {
                r: -45, g: 45, j: 135, n: -135
            }[bendtype]
            alpha1 = alpha0 - beta
            assert(-360 < alpha1 and alpha1 <= 360)
            if alpha1 <= -180: alpha1 += 360
            elif alpha1 > 180: alpha1 -= 360
            assert(-180 < alpha1 and alpha1 <= 180)
            # Normalise the cost.
            cost = abs(alpha1/180.0)
        else:
            # For a non-cycle we don't distinguish between r and j, or between g and n
            # bends, so you can only be up to +/- 90 degrees off.
            assert(-180 < alpha0 and alpha0 <= 180)
            if alpha0 <= -90: alpha0 += 180
            elif alpha0 > 90: alpha0 -= 180
            assert(-90 < alpha0 and alpha0 <= 90)
            beta = {
                r: -45, g: 45, j: -45, n: 45
            }[bendtype]
            alpha1 = alpha0 - beta
            assert(-135 < alpha1 and alpha1 <= 135)
            if alpha1 <= -90: alpha1 += 180
            elif alpha1 > 90: alpha1 -= 180
            assert(-90 < alpha1 and alpha1 <= 90)
            # Normalise the cost.
            cost = abs(alpha1/90.0)
        return cost

    def nextLocalOptimalPoint(self, i0, bendtype, remaining=0):
        """
        :param i0: a position in the chain
        :param bendtype: a bent LinkShape
        :param remaining: how many more points we must choose /after/ this one
        :return: (i1, c) being the chosen point and the cost there

        We choose a locally optimal point i1 /at or after/ position i0, at which to create
        the given bend type. Optimality means minimal cost, from the bendCost function.

        If remaining == r and there are at least r positions left after i0 in the chain,
        then we return an i1 which has at least r points left after it; if not,
        then we just return i1 = i0.
        """
        n = len(self.nodes)
        candidate = None
        bestCost = 10 # effectively infinity since costs are at most 1
        i1 = i0
        M = 2*n - 1
        if self.cycle: M += 1
        M -= remaining
        cost = bestCost
        for i in range(i0, M):
            cost = self.bendCost(bendtype, i)
            if candidate is not None and cost > bestCost:
                i1 = candidate
                cost = bestCost
                break
            # To even be considered a candidate for optimal position, the cost
            # has to be less than 0.5. Else we might start at bad and go to worse,
            # and thereby accept bad.
            if cost < 0.5 and cost < bestCost:
                candidate = i
                bestCost = cost
        else:
            if candidate is not None:
                i1 = candidate
                cost = bestCost
        return (i1, cost)

    def globalOptimalPoint(self, bendtype, beginAt=0):
        """
        :param bendtype: a bent LinkShape
        :param beginAt: a position in the chain
        :return: (i, c) being the point and the cost

        We choose a locally optimal point /at or after/ position beginAt, at which to create
        the given bend type. Optimality means minimal cost, from the bendCost function.

        If there are no points left after beginAt, we return None.
        """
        n = len(self.nodes)
        i0 = None
        cost = 10 # max cost is 1, so 10 is effectively infinity
        M = 2*n - 1
        if self.cycle: M += 1
        for i in range(beginAt, M):
            c = self.bendCost(bendtype, i)
            if c < cost:
                i0, cost = i, c
        return (i0, cost)

    def evaluateBendSeq(self, bendseq):
        """
        :param bendseq: a BendSequence object
        :return: the given bendseq object, for convenience

        We compute the best places for the prescribed bendtypes to occur and stash them in
        the bendpoints field of the bendseq object, and the cost of creating these bends in
        the cost field.

        The "places" are indices 0, 1, 2, 3, ... which refer to the first node in the chain,
        then the first edge, next node, next edge, and so on, with even numbers meaning nodes
        and odd numbers meaning edges.
        """
        queue = Deque(bendseq.bendtypes)
        i = 0
        cost = 0
        bendpoints = []
        while len(queue) > 1:
            bendtype = queue.popleft()
            i, c = self.nextLocalOptimalPoint(i, bendtype, remaining=len(queue))
            if i is not None:
                bendpoints.append(i)
                cost += c
                i += 1
        if len(queue) == 1:
            bendtype = queue.popleft()
            i, c = self.globalOptimalPoint(bendtype, beginAt=i)
            if i is not None:
                bendpoints.append(i)
                cost += c
                i += 1
        bendseq.bendpoints = bendpoints
        bendseq.cost = cost
        return bendseq

    def writeConfigSeq(self, bendseq):
        """
        :param bendseq: a BendSequence object, whose bendpoints are indices into this
                Chain's sequence of nodes AND edges -- thus even indices for nodes and
                odd indices for edges. Its corresponding bendtypes are the types of bends
                that should occur at those indices.

        :return: a "configuration sequence," which looks like
                    [ c0, c1, ..., cm-1 ]
                where m the number of edges to be configured, which is n - 1 if this is
                not a cycle, and n if it is -- n the number of nodes in the chain --
                and each ci is a list of length 1 or 2, containing Compass directions.

                When ci == [ d ], then edge i is to be configured in direction d.
                When ci == [d1, d2], then edge i is to be replaced by a bend point,
                which we go into in direction d1, and come out of in direction d2.
        """
        m = len(self.edges)
        config = []
        bends = zip(bendseq.bendpoints, bendseq.bendtypes)
        if self.cycle:
            m += 1
            assert(len(bends) == 4)
            bt0 = bends[0][1]
            # Since we always run cycles clockwise, we can infer from the first bendtype
            # what the incoming direction must be.
            dIn = cwIncomingDirForBend(bt0)
        else:
            # Not a cycle
            dIn = bendseq.incomingDirec
            assert(dIn is not None)
        ptr = 0
        direc = dIn
        for j in range(m):
            if ptr == len(bends):
                # All remaining edges get the current direc.
                config.append([direc])
            else:
                k = 2*j + 1
                bs = []
                while ptr < len(bends) and bends[ptr][0] in [k, k - 1]:
                    bs.append(bends[ptr])
                    ptr += 1
                # At this point, k is an odd number, referring to an edge in the chain,
                # direc is the incoming direction into node k - 1, and bs is a list of
                # bend points of length 0, 1, or 2, occurring at node k - 1 and/or edge k.
                # Our job is to: (a) describe what happens at edge k, namely either that
                # it be configured to run in a certain compass direction, or that it contain
                # a bend point, with certain incoming and outgoing compass directions;
                # and (b) to set direc equal to the (outgoing) direction of the edge that
                # leads into node k + 1.
                if len(bs) == 2:
                    bp0, bt0 = bs[0]
                    bp1, bt1 = bs[1]
                    dir0 = applyBendToDir(bt0, direc)
                    dir1 = applyBendToDir(bt1, dir0)
                    config.append([dir0, dir1])
                    direc = dir1
                elif len(bs) == 1:
                    bp, bt = bs[0]
                    nextDir = applyBendToDir(bt, direc)
                    assert(nextDir is not None)
                    if bp == k:
                        # Next bend should occur at this edge.
                        config.append([direc, nextDir])
                    elif bp == k - 1:
                        # Next bend should occur at the node before this edge.
                        config.append([nextDir])
                    direc = nextDir
                else:
                    # Carry on with current direction.
                    # In particular, this case handles what happens if the final
                    # bend is to occur at the final node. For in that case all we
                    # can do is carry on with the current direction, and it is
                    # up to the anchorEdgeRight to make the final bend happen.
                    config.append([direc])
        return config

    def possibleBendSeqs(self):
        """
        :return: list [s0, s1, ..., sk] where each si is a BendSequence object,
                 indicating a sequence of bends that this
                 chain may have, given its endpoints.

                 If "no bends" is a possibility, we return a BendSequence with empty
                 list of bend types.
        """
        if len(self.nodes) == 1: return None
        assert(len(self.nodes) >= 2)
        seqs = []
        if self.cycle:
            seqs = [BendSequence(LinkShape.cwBendsFrom(bt)) for bt in LinkShape.bent]
        else:
            # Get incoming and outgoing directions:
            A, Z = self.anchorNodeLeft, self.anchorNodeRight
            b, y = self.nodes[0], self.nodes[-1]
            dIn = self.graph.nodeConf.getDirec(A, b)
            dOut = self.graph.nodeConf.getDirec(y, Z)
            # If edge (A, b) or (y, Z) is not configured, we look up possible directions.
            dIns = [dIn] if dIn is not None else Compass.possibleCardinalDirections(A, b)
            dOuts = [dOut] if dOut is not None else Compass.possibleCardinalDirections(y, Z)
            # Now compute the sequences.
            for d0 in dIns:
                for d1 in dOuts:
                    seqs.extend([
                        BendSequence(bs, dIn=d0, dOut=d1)
                        for bs in lookupMinimalBendSeqs(A, d0, Z, d1)
                    ])
        return seqs

    def takeShapeBasedConfiguration(self):
        """
        Give this chain an orthogonal configuration best fitting its present
        geometric shape, i.e. putting the bend points in the most natural places,
        including the possibility that they go where edges are (meaning a new bend
        point is created).

        :return: boolean saying if anything changed
        """
        # For a chain of one node, there is nothing to do.
        if len(self.nodes) == 1: return False
        # Else there is at least one internal edge, and we assume that none of the
        # internal edges is yet configured. Therefore we /will/ be making
        # changes -- even if not creating any bent links (straight chains still need
        # to be configured).
        changes = True
        seqs = self.possibleBendSeqs()
        assert(len(seqs) > 0)
        for bs in seqs: self.evaluateBendSeq(bs)
        seqs.sort(key=lambda bs: bs.cost)
        bestSeq = seqs[0]
        configseq = self.writeConfigSeq(bestSeq)
        # Now create the configuration.
        G = self.graph
        for j, conf in enumerate(configseq):
            # Get the edge and the nodes u, v that come before and after
            # it in the chain, respectively.
            k = 2*j + 1
            edge = self.getEdge(k)
            u, v = self.getNode(k - 1), self.getNode(k + 1)
            if len(conf) == 1:
                # In this case the edge is to be aligned in a compass direction.
                G.nodeConf.setDirec(u, v, conf[0])
            else:
                assert(len(conf) == 2)
                # In this case we are to create a bend point.
                # First sever the edge and free any config on the nodes.
                G.severEdge(edge)
                G.nodeConf.free(u, v)
                # Create a bend point, giving it a reasonable initial location
                # midway between u and v. It is not yet aligned with either of
                # them; that will wait until we project onto the new constraints.
                x, y = (u.x + v.x) / 2.0, (u.y + v.y) / 2.0
                bp = G.createAndAddNewBendNode(x, y)
                bp.setIDAsLabel()
                # Configure
                G.nodeConf.setDirec(u, bp, conf[0])
                G.nodeConf.setDirec(bp, v, conf[1])
                # Add edges.
                G.createAndAddNewEdge(u.ID, bp.ID)
                G.createAndAddNewEdge(bp.ID, v.ID)
                # Save a record.
                self.aestheticBends.append(AestheticBend(bp, u, v))
        return changes

def test1():
    class DummyNode:
        def __init__(self, x, y):
            self.x = x
            self.y = y
    A = DummyNode(-30, 0)
    Z = DummyNode(0, 0)
    d0 = Compass.WEST
    d1 = Compass.NORTH
    print lookupMinimalBendSeqs(A, d0, Z, d1)

if __name__ == '__main__':
    test1()
