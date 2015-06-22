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

"""
This module is for removing edge overlaps in routed orthogonal graphs.

To say that a graph G is "routed" means that each edge in G should
have its self.route attribute defined already.

To say that G is orthogonal means that each segment of these edge routes
is axis-aligned, i.e. precisely one coordinate changes between consecutive
route points.
"""

from graphs import *
from ortho import Compass

def partition(lst, key, tol=None):
    """
    :param lst: a list
    :param key: a callable
    :param tol: a tolerance -- optional
    :return: a list of lists, representing a partition of lst
             according to the value returned by key on each element.

    If you pass a tolerance value then we use the tolerantPartition
    function instead.
    """
    if tol is not None:
        return tolerantPartition(lst, key, tol=tol)
    parts = {}
    keys = []
    for item in lst:
        k = key(item)
        if k not in keys:
            keys.append(k)
        part = parts.get(k, [])
        part.append(item)
        parts[k] = part
    partList = [
        parts[k] for k in keys
    ]
    return partList

def tolerantPartition(lst, key, tol=0.01):
    """
    Like the basic parition function, but allows key values to
    be within a given tolerance of one another.
    (Key values should be floats!)
    """
    # If lst is empty there's nothing to be done.
    if len(lst) == 0: return []
    # Else begin by sorting the list.
    lst.sort(key=key)
    # Prepare return value.
    partList = []
    # Initialise the first part with the first item of the list.
    firstItem = lst[0]
    part = [firstItem]
    # Initialise the average key value for the part.
    avg = key(firstItem)
    n = 1  # n records how many keys are in the avg
    for item in lst[1:]:
        k = key(item)
        dk = k - avg
        if abs(dk) < tol:
            # Difference of new key with current avg is within
            # tolerance. Add item to current part and update avg.
            part.append(item)
            avg = (n*avg + k)/(n+1.0)
            n += 1
        else:
            # Difference is too much. Record the current part
            # and begin a new one for this item.
            partList.append(part)
            part = [item]
            avg, n = k, 1
    # At this point there is always a nonempty part that has not
    # yet been appended to the part list.
    partList.append(part)
    return partList

class EdgeSegment:

    def __init__(self, node1, node2):
        """
        Pass the two Nodes that are the endpoints of the segment.
        We will define:
            self.kind: "H" or "V"
            self.u1, self.u2: the lower and upper bounds on the interval
                              spanned by the segment, x if H, y if V.
            self.n1, self.n2: the two endpt nodes, n1 where the segment
                              opens, n2 where it closes
            self.v: the y coord of the segment if H, the x coord if V

        Finally, Segments are not per se directed; however, it is convenient to note
        the direction implied by the order in which the arguments were passed
        to this __init__ function, and we record this in self.origDir.
        """
        x1, y1 = node1.x, node1.y
        x2, y2 = node2.x, node2.y
        dx, dy = x2 - x1, y2 - y1
        if abs(dy) <= abs(dx):
            self.kind = "H"
            self.v = y1
            if x1 < x2:
                self.u1, self.u2 = x1, x2
                self.n1, self.n2 = node1, node2
                self.origDir = Compass.EAST
            else:
                self.u1, self.u2 = x2, x1
                self.n1, self.n2 = node2, node1
                self.origDir = Compass.WEST
        else:
            self.kind = "V"
            self.v = x1
            if y1 < y2:
                self.u1, self.u2 = y1, y2
                self.n1, self.n2 = node1, node2
                self.origDir = Compass.SOUTH
            else:
                self.u1, self.u2 = y2, y1
                self.n1, self.n2 = node2, node1
                self.origDir = Compass.NORTH
        # A reference to the owning edge may be stashed here, if you wish:
        self.edge = None

    def setNewClosingNode(self, node):
        self.n2 = node
        self.u2 = node.x if self.kind == "H" else node.y

    def __repr__(self):
        k = self.kind
        v = 'y' if k == "H" else 'x'
        u = 'x' if k == "H" else 'y'
        s = '"%s"-seg at %s = %.2f over %s-interval [%.2f, %.2f]'%(
            k, v, self.v, u, self.u1, self.u2
        )
        s += ' between nodes %s and %s' % (
            self.n1.ID, self.n2.ID
        )
        return s

    def getEvents(self):
        """
        Return the two Event objects, in order, representing the
        opening and closing of this segment.
        """
        evts = [Event(self, n) for n in [self.n1, self.n2]]
        evts[0].companion = evts[1]
        evts[1].companion = evts[0]
        return evts

class Event:

    def __init__(self, seg, endpt):
        self.seg = seg
        self.endpt = endpt
        # const-coord:
        self.v = endpt.y if seg.kind == "H" else endpt.x
        # var-coord:
        self.u = endpt.x if seg.kind == "H" else endpt.y
        # kind
        if self.u == seg.u1:
            self.kind = "open"
        elif self.u == seg.u2:
            self.kind = "close"
        else:
            raise Exception("endpt does not match segment")
        self.companion = None  # for holding a ref to corresp. event

def buildSegments(G):
    """
    Before calling this method it is necessary that the routePoints
    of each edge in G be built.

    :param G: a routed orthogonal graph
    :return: a list of EdgeSegment objects, one for each segment in G
    """
    segs = []
    for edge in G.edges.values():
        pts = [edge.src] + edge.routePoints + [edge.tgt]
        n = len(pts)
        for i in range(n-1):
            p, q = pts[i:i+2]
            seg = EdgeSegment(p, q)
            seg.edge = edge
            segs.append(seg)
    return segs

def computeNodeGroups(segs):
    """
    :param segs: a list of Segments (should be all H- or all V-segs)
    :return: a list of lists of Nodes ("groups"), being endpts that
             participate in a common sequence of overlapping segments,
             listed in order of increasing variable-coordinate
    """
    gps = []
    # Partition segs by const-coord
    #parts = partition(segs, lambda s: s.v, tol=0.01)
    parts = partition(segs, lambda s: s.v, tol=0.5)
    # Compute groups for each part.
    for part in parts:
        # Build list of events.
        evts = []
        for seg in part:
            evts.extend(seg.getEvents())
        # Sort events by variable-coord
        evts.sort(key=lambda e: e.u)
        # Initialise
        gp = []
        openSegs = set([])
        for e in evts:
            endpt = e.endpt
            # Do not append multiple copies:
            if len(gp) == 0 or gp[-1] is not endpt:
                gp.append(endpt)
            # Keep track of which segs are open.
            if e.kind == 'open':
                openSegs.add(e.seg)
            elif e.kind == 'close':
                openSegs.remove(e.seg)
                # If no open segs remain, the group is complete.
                if len(openSegs) == 0:
                    gps.append(gp)
                    gp = []
    return gps

def removeEdgeOverlaps(G):
    """
    Remove the edge overlaps in a routed orthogonal graph G.
    :param G: The graph whose overlaps are to be removed.
    :return: A new graph object Q. It has all the nodes of G, plus
             a node for each bend point in the routes of G. Its edges
             cover the routes of G, in a set-theoretic sense, but none
             of its edges overlap.
    """
    DEBUG = False
    # Get the new graph Q started by giving it a copy of each node in G,
    # and a node for each bend point in G.
    Q = Graph()
    for node in G.nodes.values():
        cnode = node.excisedCopy()
        Q.addNode(cnode)
    bps = G.buildBendPoints()
    for node in bps.values():
        Q.addNode(node)
    # Build the EdgeSegment objects for the edges in G.
    segs = buildSegments(G)
    # Testing:
    if DEBUG:
        for seg in segs:
            print repr(seg)
    # Separate into H-segs and V-segs
    hsegs = filter(lambda s: s.kind == "H", segs)
    vsegs = filter(lambda s: s.kind == "V", segs)
    # Compute node groups
    hgps = computeNodeGroups(hsegs)
    vgps = computeNodeGroups(vsegs)
    # Testing:
    if DEBUG:
        print 'H-groups:'
        for g in hgps:
            print [n.ID for n in g]
        print 'V-groups:'
        for g in vgps:
            print [n.ID for n in g]
    # Create edges for Q based on the node groups.
    gps = []
    gps.extend(hgps)
    gps.extend(vgps)
    for gp in gps:
        for i in range(len(gp) - 1):
            s, t = gp[i:i+2]
            e = Edge(s.ID, t.ID)
            Q.addEdge(e)
    #if DEBUG:
    #    Q.setIDsAsLabels()
    return Q

def computeCrossings(segs, idDisp):
    """
    ...
    :param segs:
    :return:
    """
    crossingNodes = {}
    # List all events for all segs.
    evts = []
    for seg in segs:
        evts.extend(seg.getEvents())
    # Sort and partition by x.
    xkey = lambda e: e.endpt.x
    evts.sort(key=xkey)
    #xparts = tolerantPartition(evts, xkey)
    #xparts = partition(evts, xkey)
    xparts = partition(evts, xkey, tol=0.8)

    #xparts = partition(evts, xkey, tol=1.0)
    # Scan through by x-coord.

    # How to choose setting for SORT_BY_EXACT_Y and TOLERANCE:
    # Suppose vertical segment A has its south end at (0, 0), and horizontal
    # segment B has its east end at (0, -0.00000000001). This means that
    # /technically/ A and B intersect. However (http://xkcd.com/1475/) you
    # probably don't actually want to treat this as an intersection. In that
    # case, set SORT_BY_EXACT_Y to False. If it doesn't help, consider making
    # the TOLERANCE bigger.
    #
    # In this example, setting SORT_BY_EXACT_Y to False means that, when the
    # list of active events is sorted, the "close" event for segment A will come
    # /before/ the "sustain" event for segment B, instead of the other way around,
    # as dictated by their exact y-coordinates. This way we will /not/ detect an
    # intersection between A and B.
    SORT_BY_EXACT_Y = False
    #TOLERANCE = 0.0001
    #TOLERANCE = 0.1
    TOLERANCE = 1.0

    def activeCmp(e,f):
        """
        Sort primarily by increasing y-coord.
        Secondarily, "close" events come before "sustain," and those
        before "open" events.
        """
        EPSILON = 0 if SORT_BY_EXACT_Y else TOLERANCE
        if f.endpt.y - e.endpt.y > EPSILON:
            return -1
        elif e.endpt.y - f.endpt.y > EPSILON:
            return 1
        else:
            kindNums = {
                'close': 0,
                'sustain': 1,
                'open': 2
            }
            ek = kindNums[e.kind]
            fk = kindNums[f.kind]
            return ek - fk

    openH = set([])
    for part in xparts:
        openV = None
        active = []
        active.extend(openH)
        active.extend(part)
        active.sort(cmp = activeCmp)
        for evt in active:
            if evt.kind == 'sustain' and openV is not None:
                # Create new crossing node.
                cn = Node()
                cn.ID = idDisp.takeNext()
                cn.x, cn.y = openV.v, evt.v
                cn.w, cn.h = 8, 8
                cn.shape = 'circle'
                cn.fill = '#FF0000'
                cn.setIDAsLabel()
                crossingNodes[cn.ID] = cn
                # Set cn as new endpt of the two segs.
                evt.seg.setNewClosingNode(cn)
                openV.seg.setNewClosingNode(cn)
                # New H-seg
                hseg = EdgeSegment(cn, evt.companion.endpt)
                segs.append(hseg)
                evt.companion.seg = hseg
                # new V-seg
                vseg = EdgeSegment(cn, openV.companion.endpt)
                segs.append(vseg)
                openV.companion.seg = vseg
                # update the two open events
                # H:
                evt.seg = hseg
                evt.endpt = cn
                evt.u = cn.x
                # and simply allow evt to stay in openH for the next part
                # V:
                openV.seg = vseg
                openV.endpt = cn
                openV.u = cn.y
                # and simply allow openV to remain the open vertical event
            elif evt.kind == 'open':
                if evt.seg.kind == "H":
                    evt.kind = 'sustain'
                    openH.add(evt)
                elif evt.seg.kind == "V":
                    openV = evt
            elif evt.kind == 'close':
                if evt.seg.kind == "H":
                    openH.remove(evt.companion)
                elif evt.seg.kind == "V":
                    openV = None
    return crossingNodes

def removeEdgeCrossings(Q, withConstraints=False):
    """
    :param Q: a routed orthogonal graph /with no edge overlaps/
    :param withConstraints: say whether you want constraints to be
                            generated for each edge
    :return: a planarisation P of Q, obtained by inserting dummy
             nodes at all edge crossings in Q
    """
    DEBUG = False
    # Get the new graph P started by giving it a copy of each node in Q.
    P = Graph()
    for node in Q.nodes.values():
        cnode = node.excisedCopy()
        P.addNode(cnode)
    # Build a EdgeSegment object for each edge in Q.
    segs = buildSegments(Q)
    # Compute the crossing pairs.
    crossingNodes = computeCrossings(segs, Q.getIDdispenser())
    # Testing:
    if DEBUG:
        print 'Crossing nodes:'
        for ID in crossingNodes:
            cn = crossingNodes[ID]
            print '%s: %s' % (ID, cn)
        for seg in segs:
            print seg
    for cn in crossingNodes.values():
        P.addNode(cn)
    for seg in segs:
        s, t = seg.n1, seg.n2
        e = Edge(s.ID, t.ID)
        P.addEdge(e)
    if withConstraints:
        for seg in segs:
            d = EdgeSegment(seg.n1, seg.n2).origDir
            P.nodeConf.setDirec(seg.n1, seg.n2, d)
    return P
