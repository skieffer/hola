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

from logging import LogLevel
from collections import deque as Deque
from ortho import *
from graphs import newBendNode
from graphs import Edge
from graphs import Node
from graphs import Graph
from graphs import RoutingRig
from constraints import SepCo, ProjSeq
from chainconfig import *
from water import WaterDivide
from adaptagrams.adaptagrams import XDIM
from adaptagrams.adaptagrams import YDIM
from adaptagrams.adaptagrams import OrthogonalRouting


class FaceSet:
    """
    Holds all the Face objects for a given Graph, and provides methods to
    use and manage them.
    """

    def __init__(self, G, logger, config):
        """
        :param G: a Graph

        We compute all the faces of the graph, storing each as a Face object,
        we determine which is the external face,
        we index the faces by the nodes they contain.
        """
        self.graph = G
        self.logger = logger
        self.config = config
        self.faces = self.computeFaces()
        self.externalFace = self.identifyExternalFace()
        # Index the faces by the nodes they contain.
        self.nodeIDtoFaces = {}
        for F in self.faces:
            for node in F.nodeseq:
                fs = self.nodeIDtoFaces.get(node.ID, [])
                fs.append(F)
                self.nodeIDtoFaces[node.ID] = fs

    def getAllTreePlacements(self):
        tps = []
        for F in self.faces:
            tps.extend(F.getAllTreePlacements())
        return tps

    def applyGeometryToTrees(self, iel=0):
        for F in self.faces:
            F.applyGeometryToTrees(iel=iel)

    def numTreesWithGrowthDirec(self, dg, scaleBySize=False):
        return sum(F.numTreesWithGrowthDirec(dg, scaleBySize=scaleBySize)
                   for F in self.faces)

    def insertTreesIntoGraph(self, H):
        """
        :param H: a Graph
        :return: nothing

        Add the individual nodes of the trees to the graph H.
        """
        for F in self.faces:
            F.insertTreesIntoGraph(H)

    def generateAllTreeSepcos(self):
        """
        :return: list of the tree SepCos for each face
        """
        sepcos = []
        for F in self.faces:
            sepcos.extend(F.generateAllTreeSepcos())
        return sepcos

    def listAllPossibleTreePlacements(self, tree):
        """
        :param tree: a Tree to be placed
        :return: a list of TreePlacement objects giving all the possible placements
                 of the tree at its root node, into any of the faces that contain that
                 node
        """
        # We actually need the node in self.graph whose ID equals that of the tree's
        # root, not the tree's root itself.
        root = self.graph.nodes[tree.root.ID]
        # Get list of available faces.
        facesAvail = self.nodeIDtoFaces[root.ID]
        # Get list of all possible TreePlacements.
        tps = []
        for F in facesAvail:
            tps.extend(F.allPossibleTreePlacements(tree, root))
        return tps

    def identifyExternalFace(self):
        """
        :return: the external face

        We identify the external face and return it. We also set each Face object's
        self.external boolean.

        NB: We assume G is 4-planar orthogonal; though it is not hard
        to do this in full generality, for now we put that off.

        To do it in full generality: instead of grabbing a node of maximal x-coord,
        get one of minimal x-coord. Then among the neighbours of this node you just
        want the one whose atan2 w.r.t. the first node is maximal. Simple.
        (The range of atan2 is (-pi, pi] over the domain from infinitessimally north
         of west, clockwise to exactly west.)
        """
        G = self.graph
        extface = None
        # Begin by grabbing any node u in G of maximal x-coord.
        u = sorted(G.nodes.values(), key=lambda v: v.x)[-1]
        # The node u cannot have any neighbour to the east. Therefore the node v following
        # u in the traversal of the exterior face must be that which comes first in
        # anticlockwise order among u's neighbours, starting with north. Since the degree
        # of u is at least two, and it can have at most one south neighbour, it must have
        # at least one neighbour among north and west. Therefore v is u's north neighbour
        # if it has one, else u's west neighbour.
        for edge in u.edges.values():
            v = edge.otherEnd(u)
            d = Compass.cardinalDirection(u, v)
            if d == Compass.NORTH: break
            elif d == Compass.WEST: w = v
        else:
            v = w
        sig = (u.ID, v.ID)
        # The external face is then the unique one containing the "signature" sig as a
        # subsequence in its sequence of node IDs (considered cyclically).
        for F in self.faces:
            IDs = [node.ID for node in F.nodeseq]
            try:
                i0 = IDs.index(sig[0])
            except:
                continue
            else:
                i1 = (i0 + 1) % F.n
                if IDs[i1] == sig[1]:
                    F.external = True
                    extface = F
                    break
        return extface

    def computeFaces(self):
        """
        :return: list of Face objects, one for each face in the graph G

        Face traversal code adapted from
            http://mathoverflow.net/questions/23811/reporting-all-faces-in-a-planar-graph
        which in turn is an adaptation of routine from SAGE:
            http://www.sagemath.org/doc/reference/graphs/sage/graphs/generic_graph.html#sage.graphs.generic_graph.GenericGraph.faces

        As per discussion at the Mathoverflow page, it is a pretty simple idea:
           "
           I'll assume the graph is connected, and that you have the clockwise or
           counterclockwise ordering of the edges around each vertex. Then it's easy,
           given a directed edge e, to walk around the face whose counterclockwise boundary
           contains e. So make a list of all directed edges (i. e., two copies of each
           undirected edge). Pick one directed edge, walk counterclockwise around its face,
           and cross off all the directed edges you traverse. That's one face. Pick a
           directed edge you haven't crossed off yet and walk around its face the same way.
           Keep doing that until you've crossed off all of the edges. (Note that the
           "counterclockwise" boundary of the exterior unbounded face actually goes clockwise
           around the outside of the graph.)
           "
        """
        G = self.graph
        # Build a /list/ of directed edges, represented as ordered pairs of node
        # IDs, with two directed edges per edge of the given graph G, one in
        # each direction.
        #
        # Used to use a set here instead of a list, but it is generally agreed upon:
        #     http://stackoverflow.com/questions/3406341/iteration-order-of-sets-in-python
        # that counting on a consistent iteration order over a set is a bad idea.
        # We sacrifice a modicum of speed in an effort to eradicate all randomness.
        edgelist = []
        for dsc in sorted(G.edges.keys()):
            edge = G.edges[dsc]
            s, t = edge.srcID, edge.tgtID
            edgelist.extend([(s, t), (t, s)])
        # Build embedding, a dictionary that looks like:
        #     v1:[v2,v3], v2:[v1], v3:[v1]
        # giving the clockwise ordering of neighbors at each vertex.
        # (Here the vi are IDs of nodes.)
        embedding = {}
        for ID in sorted(G.nodes.keys()):
            node = G.nodes[ID]
            nbrs = node.getNbrsCWCyclic()
            embedding[node.ID] = [n.ID for n in nbrs]
        faces = []
        # Initialise the first path.
        path  = []
        path.append(edgelist.pop())
        # Trace faces
        while (len(edgelist) > 0):
            # Get nbrs of target node of final edge e' in path.
            nbrs = embedding[path[-1][-1]]
            # The source node of edge e' is among those neighbours.
            # Get its index in the list of neighbours, and decrement it
            # (mod number of nbrs), in order to traverse the face in
            # clockwise order.
            next_index = (nbrs.index(path[-1][-2]) - 1) % (len(nbrs))
            next_node = nbrs[next_index]
            # Form the next edge in the path.
            tup = (path[-1][-1],next_node)
            # If it equals the first edge in the path, we have finished
            # traversing this face.
            if tup == path[0]:
                F = Face(G, self.logger, self.config)
                F.initWithEdgeSeq(path)
                F.ID = len(faces)
                faces.append(F)
                path = []
                path.append(edgelist.pop())
            # Else add it to the current path and continue.
            else:
                path.append(tup)
                edgelist.remove(tup)
        if (len(path) != 0):
            F = Face(G, self.logger, self.config)
            F.initWithEdgeSeq(path)
            F.ID = len(faces)
            faces.append(F)
        return faces


class TreePlacement:

    def __init__(self, tree, face, trunkroot, dp, dg, flip=False):
        # the Tree object:
        self.tree = tree
        # the Face into which the Tree is to be placed:
        self.face = face
        # a Node belonging to the Face which is to be the root of the Tree:
        self.node = trunkroot
        # the Compass direc in which the Tree is to be placed:
        self.placementDirec = dp
        # the Compass direc in which the Tree is to grow:
        self.growthDirec = dg
        # bool saying if the Tree is to be flipped:
        self.flip = flip
        # a float giving the cost of the placement, by some metric:
        self.cost = 0
        # place for a list of ProjSeqs, giving ways to expand the face to make
        # room for the Tree:
        self.expansionOptions = None
        # place to hold a Node representing this tree's bounding box:
        self.boxNode = None

    def computeAlignedSets(self):
        """
        :return: dictionary {'h': Hsets, 'v': Vsets}, where Hsets is a list of lists,
                 each memeber of which is a list of nodes that are horizontally aligned;
                 likewise with Vsets.

                 This placement's root is substituted for this tree's root node.
        """
        return self.tree.computeAlignedSets(self.node)

    def buildBufferNodesAndPCs(self, iel, doBuildBufferNodes=True):
        """
        This is for building the simpler set of buffer nodes. See method
        by same name in Tree class.
        """
        def makeNode(x, X, y, Y):
            u = Node()
            u.w, u.h = X - x, Y - y
            u.x, u.y = x + u.w/2.0, y + u.h/2.0
            u.ID = self.face.graph.getNextID()
            u.fill = '#C0804080'
            u.setIDAsLabel()
            return u
        bns = []
        extraNodes = {}
        if doBuildBufferNodes:
            pad = iel/4.0
            x, X, y, Y = self.boxNode.boundingBoxxXyY()
            dg = self.growthDirec
            if dg in Compass.horizontal:
                bns.append(makeNode(x, X, y, y + pad))
                bns.append(makeNode(x, X, Y - pad, Y))
                if dg == Compass.EAST:
                    bns.append(makeNode(X - pad, X, y + pad, Y - pad))
                else:
                    assert(dg == Compass.WEST)
                    bns.append(makeNode(x, x + pad, y + pad, Y - pad))
            else:
                assert(dg in Compass.vertical)
                bns.append(makeNode(x, x + pad, y, Y))
                bns.append(makeNode(X - pad, X, y, Y))
                if dg == Compass.SOUTH:
                    bns.append(makeNode(x + pad, X - pad, Y - pad, Y))
                else:
                    assert(dg == Compass.NORTH)
                    bns.append(makeNode(x + pad, X - pad, y, y + pad))
            extraNodes = {u.ID: u for u in bns}
        frc = self.tree.writeFixedRelCo(self.node, extraNodes=extraNodes)
        pcs = [frc]
        return (bns, pcs)

    def insertTreeIntoGraph(self, H, iel):
        """
        :param H: a Graph
        :param iel: ideal edge length for the graph
        :return: nothing

        Insert the nodes of the tree into the graph H.
        Remove the treebox from H if it was present.
        Build buffer nodes for the tree, and add them to H.
        Write constraints to maintain the shape of the tree, and to
        keep the buffer nodes beside the proper nodes.
        """
        self.bufferNodes = []
        # Insert the tree nodes.
        self.tree.insertIntoGraph(H)
        # Remove tree box.
        B = self.boxNode
        if H.nodes.has_key(B.ID):
            H.severNodes({B.ID:B})
        # Build buffer nodes and write PCs.
        bns, pcs = self.tree.buildBufferNodesAndPCs(iel, self.placementDirec)
        self.bufferNodes = bns
        for bn in bns:
            H.addNode(bn)
        H.addPCHolders([self])

    def rotatePCsCW(self, n, iel):
        """
        Part of PCHolder interface.
        """
        self.tree.growthDir = Compass.rotateCW(2*n, self.tree.growthDir)
        self.tree.buildBufferNodesAndPCs(iel, self.placementDirec)

    def getCCs(self, rs=None, ix=None):
        """
        Part of PCHolder interface. Generate cola constraints, and return them.
        :return: list of cola constraints
        """
        return self.tree.getCCs(rs=rs, ix=ix)

    def applyGeometryToTree(self, iel=0):
        """
        Rotate, flip, and translate the tree as necessary to match this placement.
        :return: nothing
        """
        dg = self.growthDirec
        self.tree.rotateTo(dg)
        if self.flip:
            self.tree.flip(dg)
        vect = self.node.centre()
        if self.placementDirec in Compass.cwCards:
            self.tree.translate(vect, dg)
        else:
            x0, y0, w0, h0, u0, v0 = self.tree.getBoundingBoxXYWHWithoutRoot(
                                        growthDir=dg)
            w1, h1, u1, v1 = self.treeBoxWithRootVector(iel=iel)
            dw = -iel/8.0 if dg in Compass.increasing else iel/8.0
            if dg in Compass.vertical:
                v1 += dw
            else:
                u1 += dw
            x, y = vect
            vect = (x + u1 - u0, y + v1 - v0)
            self.tree.translate(vect, dg)
        # Do this last, after we have rotated the tree.
        self.tree.growthDir = dg

    def __repr__(self):
        s = ''
        s += 'Tree placement: at node %d, placed %s' % (
            self.node.ID, Compass.abbrev[self.placementDirec]
        )
        if self.placementDirec in Compass.cwCards:
            if self.flip:
                s += ', flipped'
        else:
            s += ', growing %s' % Compass.abbrev[self.growthDirec]
        if self.expansionOptions is not None:
            s += ': cost = %.3f' % self.cost
        return s

    def getNumPotentialNbrs(self):
        """
        :return: number of "potential neighbours" of this tree if placed according to
                 this tree placement -- this is equal to the number of other root nodes
                 on the relevant Sides to which the root node of this placement belongs
        """
        sides = self.face.getRelevantSidesForPlacement(self)
        # Total number of root nodes present:
        nr = sum([S.getNumRootNodes() for S in sides])
        # Own root node is counted once per side, so subtract that.
        return nr - len(sides)

    def computeBoundarySegments(self, iel):
        """
        :param iel: ideal edge length
        :return: list of BoundarySegment objects, representing the part of the
                 interior face boundary made up by the polygon that is the union of
                 the treenode with the rootnode
        """
        pts = self.polygonPoints(iel)
        # We will want a function that tells us how the path bends at a given point.
        def bendDir(i):
            "return -1, 0, 1 if the bend is left, straight, right"
            p = pts[-1:] + pts + pts[:1]
            a, b, c = p[i:i+3]
            # If either of the neighbours a, c of b lies at the same point as b
            # (this happens in the case of a bridge edge), then we simply return 0,
            # indicating that a bend does not happen at b.
            ax, ay = a
            bx, by = b
            cx, cy = c
            if ax == bx and ay == by: return 0
            if cx == bx and cy == by: return 0
            # Otherwise b differs from both of its neighbours, and we proceed to
            # check whether it is a bend point.
            d = Compass.cardinalDirection(a, b)
            e = Compass.cardinalDirection(b, c)
            di = Compass.cwCards.index(d)
            ei = Compass.cwCards.index(e)
            r = (ei+4-di)%4
            r = -1 if r == 3 else r
            assert(r in [-1, 0, 1])
            return r
        # The points of the polygon include the points where edges connect the root
        # node to its neighbours in the face, and these neighbours may be four or two
        # in number, depending on whether the root node is or is not a pinch node,
        # respectively.
        #
        # Also there may be three neighbours if one of the edges is a bridge edge,
        # or two neighbours with two bridge edges.
        #
        # In all cases there are either two or four edgepoints, and we have sorted the
        # list of points so that the very last point is the final edgepoint.
        # Moreover the points participating in the /interior/ boundary of the face
        # are those between pairs of edgepoints, taken in order.
        #
        # I.e. if there are precisely two edge points then we want the segs between these.
        # Otherwise there are four, and we want the segs between the 1st and 2nd, as
        # well as the segs between the 3rd and the 4th.
        interiorPointLists = []
        pairs = self.face.nbrPairs[self.node.ID]
        if len(pairs) == 1:
            i0 = 0
            while bendDir(i0) != 0: i0 += 1
            interiorPointLists.append(pts[i0:])
        else:
            assert(len(pairs) == 2)
            i0 = 0
            while bendDir(i0) != 0: i0 += 1
            i1 = i0 + 1
            while bendDir(i1) != 0: i1 += 1
            i2 = i1 + 1
            while bendDir(i2) != 0: i2 += 1
            interiorPointLists.append(pts[i0:i1+1])
            interiorPointLists.append(pts[i2:])
        # It is only at concave bend points that segment endpoints represent
        # a crossing btw interior and exterior of the face.
        # Since the segments are oriented so that the left-hand side is the inside
        # of the polygon, the concave bends are precisely where we turn /right/
        # when going from one segment to the next.
        segs = []
        for pts in interiorPointLists:
            N = len(pts)
            for i in range(N-1):
                p0, p1 = pts[i:i+2]
                x0 = (i != 0 and bendDir(i) == 1)
                x1 = (i+1 != N-1 and bendDir(i+1) == 1)
                segs.append(BoundarySegment(p0, x0, p1, x1, [self]))
        return segs

    def polygonPoints(self, iel):
        """
        :param iel: ideal edge length
        :return: list of points (x, y) of the polygon that is the union of
                 the treenode with the rootnode

        If the placement direction is ordinal, then we make a slight extension to
        the treebox so that the two boxes do not meet at a mere point.

        The polygon points proceed anticlockwise around the rootnode+treenode combination,
        so that the right-hand side of each segment is the interior of the face,
        IF the segment is inside the face. The left-hand side of each segment is the
        interior of the polygon.

        We add an extra point where each edge of the face connects to the
        root node. These will be the only points whose two neighbouring segments
        run in the same dimension as each other. We call these points "edgepoints."

        The points are ordered so that the segments belonging to the /interior/ boundary
        of the face are those from the first to second edgepoints inclusive, as well as
        those from the third to fourth edgepoints inclusive in the case that the latter
        exist, i.e. in the case that the rootnode is a pinch node.
        """
        ISTHMUS_WIDTH = 4
        EPSILON = 0.01
        rx, ry = self.node.x, self.node.y
        w, h, u, v = self.treeBoxWithRootVector(iel=iel)
        dp, dg = self.placementDirec, self.growthDirec
        # Add "isthmus" for ordinal placement directions.
        if dp in Compass.cwOrds:
            iw = ISTHMUS_WIDTH
            if dg in Compass.horizontal:
                w += iw
                u += -iw/2.0 if dg in Compass.increasing else iw/2.0
            else:
                h += iw
                v += -iw/2.0 if dg in Compass.increasing else iw/2.0
        # Make dummy node for treebox
        treebox = Node()
        treebox.x, treebox.y = rx + u, ry + v
        treebox.w, treebox.h = w, h
        # Determine the side of the root on which it meets the treebox.
        if dp in Compass.cwCards:
            contactSide = dp
        else:
            contactSide = filter(lambda d: d != dg, Compass.components[dp])[0]
        # The corner of the root where we transfer to the treebox is the first
        # corner of this side, when going anticlockwise.
        transitionCorner = {
            Compass.EAST: Compass.SE, Compass.SOUTH: Compass.SW,
            Compass.WEST: Compass.NW, Compass.NORTH: Compass.NE
        }[contactSide]
        # We will visit the corners of the boxes in anticlockwise order.
        remRootDirecs = Deque(Compass.acwOrds[:])
        d0 = remRootDirecs.popleft()
        # Build list of dirs on root that we hit before hitting the tree box.
        dirs = []
        dirs.append(d0)
        # Now work anticlockwise through the ordinal directions until hitting
        # the transition corner.
        while not dirs[-1] == transitionCorner:
            dirs.append(remRootDirecs.popleft())
        # Initialise the list of points.
        pts = [self.node.getBdryCompassPt(d) for d in dirs]
        # Get first corner for tree box.
        # This is the corner 90 degrees clockwise from the transition corner.
        # (Draw the pictures and see for yourself.)
        d1 = Compass.cw90(transitionCorner)
        # Go through all the tree box points starting from here.
        aa = Compass.acwOrds + Compass.acwOrds
        treeBoxDirecs = aa[aa.index(d1): aa.index(d1)+4]
        pts.extend([treebox.getBdryCompassPt(d) for d in treeBoxDirecs])
        # And finish up with any remaining root node directions.
        pts.extend([self.node.getBdryCompassPt(d) for d in remRootDirecs])
        # Now if a side of the root node and a side of the tree box happen to be
        # in flush alignment with each other, then there will be a doubled point
        # in the points list. We now delete any such doubles (delete /both/ points),
        # thereby fusing those flush sides into one long side.
        # Due to floating point rounding errors, we detect such "identical" points
        # using a small epsilon.
        loop = pts + pts[:1]
        N, i = len(loop), 0
        pts = []
        def samePts(p, q):
            px, py = p
            qx, qy = q
            return abs(px - qx) < EPSILON and abs(py - qy) < EPSILON
        while i + 1 < N:
            p, q = loop[i:i+2]
            if samePts(p, q):
                i += 2
            else:
                pts.append(p)
                i += 1
        # Finally we insert the extra points for the edge directions.
        # Get the neighbour pairs for the root node. There may be two pairs,
        # if the root node happens to be a pinch node.
        pairs = self.face.nbrPairs[self.node.ID]
        # Now we consider the edgepoints in the direction of these neighbouring
        # nodes. Since in each neighbour pair the nodes are listed in the order of
        # clockwise traversal of the face edges, the corresponding edgepoints on
        # the perimeter of the root node are such that this node has no other
        # edgepoints between the first and the second in the anticlockwise
        # direction.
        # For example, if a pair of neighbouring nodes lie to the south, and to
        # the west, in that order, then this node must not have neighbours to
        # the east or north participating in this same face.
        # Moreover, all the directions from that of the first neighbour in a pair
        # to that of the second must be interior directions, pointing into the
        # face.
        epts = Deque()
        for pair in pairs:
            for nbr in pair:
                d = Compass.cardinalDirection(self.node, nbr)
                epts.append(self.node.getBdryCompassPt(d))
        # Insert each edge point where it belongs.
        def liesBetween(b, a, c):
            ax, ay = a
            bx, by = b
            cx, cy = c
            if abs(ax - cx) < EPSILON:
                w, wl, wh = by, min(ay, cy), max(ay, cy)
                z, z0 = bx, (ax + cx)/2.0
            elif abs(ay - cy) < EPSILON:
                w, wl, wh = bx, min(ax, cx), max(ax, cx)
                z, z0 = by, (ay + cy)/2.0
            else:
                raise Exception('pts a and c not axis aligned')
            return (abs(z - z0) < EPSILON) and (wl <= w <= wh)
        i = 0
        n = len(pts)
        while len(epts) > 0:
            ept = epts.popleft()
            i0 = i
            while True:
                a, c = pts[i], pts[(i+1)%n]
                if liesBetween(ept, a, c):
                    i1 = i + 1
                    # i1 is deliberately /not/ reduced mod n
                    # This means if i1 == n then ept is simply appended at the end of the list.
                    pts.insert(i1, ept)
                    n += 1
                    i = (i+2) % n
                    # Yes, we deliberately reduce mod the /new/ n.
                    break
                else:
                    i = (i+1) % n
                    if i == i0:
                        raise Exception('nowhere to place edgepoint')
        # Finally, ensure that the last edge point is the very last point in the list.
        # This ensures that the points of the polygon that are /interior/ to the face
        # are those between the first and second edgepoints, and those between the
        # third and fourth edgepoints, if the latter exist.
        pts = pts[i1+1:] + pts[:i1+1]
        return pts

    def somePointOppositeSegment(self, seg, iel=0, openInterval=False):
        """
        :param seg: a LineSegment
        :param iel: ideal edge length -- used to compute padding on treenode
        :param openInterval: say if we should consider the open interval of the segment
        :return: a point (x, y) or None

        We compute the interval I of the treenode in the dimension parallel to the
        segment. Let J be the interval of the segment -- open if the openInterval
        argument is True; closed otherwise. Let K be the intersection of I and J.
        If K is empty we return None. Otherwise we pick a value w in K and return
        a point having w as one of its coordinates and a centre coordinate of the
        treenode as the other coordinate.
        """
        w, h, u, v = self.treeBoxWithRootVector(iel=iel)
        rcx, rcy = self.node.x, self.node.y
        tcx, tcy = rcx + u, rcy + v
        x, X, y, Y = tcx - w/2.0, tcx + w/2.0, tcy - h/2.0, tcy + h/2.0
        I = [x, X] if seg.direc in Compass.horizontal else [y, Y]
        if openInterval:
            K = seg.openIntervalIntersection(I)
        else:
            K = seg.closedIntervalIntersection(I)
        pt = None
        if K is not None:
            w = K[0]
            pt = (w, tcy) if seg.direc in Compass.horizontal else (tcx, w)
        return pt

    def treeBoxWithRootVector(self, iel=0):
        """
        :param iel: ideal edge length for the graph. Set this to a positive value if
                    you want padding based on this length to be added to the box
                    dimensions.
        :return: (w, h, u, v) where w, h are the width and height of the bounding box
                 for the tree (minus root) for the growth direction of this TreePlacement,
                 and (u, v) is the vector from the centre of the trunk root to the centre
                 of this box, for the placement direction and flip of this TreePlacement.

                 If iel is positive, we add padding to w and h.
        """
        x0, y0, w0, h0, u0, v0 = self.tree.getBoundingBoxXYWHWithoutRoot(growthDir=self.growthDirec)
        #w0, h0, u0, v0 = self.tree.getSubRootBox(self.growthDirec)
        w1, h1 = w0, h0
        if self.placementDirec == self.growthDirec:
            # Placement direction is cardinal and equals the growth direction.
            assert(self.placementDirec in Compass.cwCards)
            # In this case we need to check whether the tree is to be flipped,
            # and if so alter the (u, v) vector accordingly.
            # Also add padding if iel is set.
            if self.growthDirec in Compass.vertical:
                u1, v1 = (-u0, v0) if self.flip else (u0, v0)
                if iel > 0:
                    w1, h1 = w1 + iel/2.0, h1 + iel/4.0
                    # Want half the width padding on each side, but all height padding on outside:
                    v1 += iel/8.0 if self.growthDirec in Compass.increasing else -iel/8.0
            else:
                assert(self.growthDirec in Compass.horizontal)
                u1, v1 = (u0, -v0) if self.flip else (u0, v0)
                if iel > 0:
                    w1, h1 = w1 + iel/4.0, h1 + iel/2.0
                    # Want half the height padding on each side, but all width padding on outside:
                    u1 += iel/8.0 if self.growthDirec in Compass.increasing else -iel/8.0
        else:
            # Placement direction is ordinal, and growth direction is one of its components.
            assert(self.placementDirec in Compass.cwOrds)
            assert(self.growthDirec in Compass.components[self.placementDirec])
            # In this case we throw away the given (u0, v0), and compute this vector
            # based solely on the dimensions of the tree and root, and on the placement direction.
            # It does not matter if the tree is to be flipped.
            # First consider padding.
            if iel > 0:
                if self.growthDirec in Compass.vertical:
                    w1, h1 = w1 + iel/2.0, h1 + iel/4.0
                else:
                    w1, h1 = w1 + iel/4.0, h1 + iel/2.0
                #w1 += iel/4.0
                #h1 += iel/4.0
            rootW, rootH = self.node.w, self.node.h
            sgnX, sgnY = Compass.vectorSigns(self.placementDirec)
            u1, v1 = sgnX * (rootW + w1)/2.0, sgnY * (rootH + h1)/2.0
        return (w1, h1, u1, v1)


    def __lt__(self, other):
        # Cardinal placement direc beats ordinal.
        sc = self.placementDirec in Compass.cwCards
        oc = other.placementDirec in Compass.cwCards
        if sc and not oc:
            return True
        elif oc and not sc:
            return False
        else:
            # Otherwise external face beats internal.
            se = self.face.external
            oe = other.face.external
            if se and not oe:
                return True
            elif oe and not se:
                return False
            else:
                # Otherwise compare cost.
                if self.cost < other.cost:
                    return True
                else:
                    return False

    def evaluateExpansionOptions(self, iel):
        """
        In general we may need to expand the face to make room for the tree.
        This method asks self.face to work out all the ways to expand, try
        them out, and sort them by how much they increase the stress of the graph.
        We store the options in self.expansionOptions as a list of ProjSeq objects.

        :param iel: we need to know the ideal edge length for the graph in order
                    to evaluate stress; also to put padding around tree boxes
        :return: nothing
        """
        self.expansionOptions = self.face.evaluateExpansionOptions(self, iel)
        self.cost = self.expansionOptions[0].getTotalStressChange()

    def getBestProjSeq(self, iel):
        """
        :param iel: ideal edge length for the graph
        :return: the best expansion option, as a ProjSeq object
        """
        if self.expansionOptions is None:
            self.evaluateExpansionOptions(iel)
        return self.expansionOptions[0]

    def evaluateCost(self, iel):
        """
        :param iel: ideal edge length for the graph
        :return: cost

        We evaluate the expansion options if that has not been done already,
        store the results, and return the cost of the best one.
        """
        if self.expansionOptions is None:
            self.evaluateExpansionOptions(iel)
        return self.cost

    def estimateCost(self, iel, logger):
        # For this purpose the primary dimension does not matter, so
        # just pass XDIM.
        wd = WaterDivide(self, XDIM, iel, logger)
        self.cost = wd.estimateCost()
        return self.cost

    def insertTreeNode(self, iel):
        """
        To be used after the face has been expanded to make room for the tree.
        This method adds a large node to the graph, representing the bounding
        box of the tree. The treenode is constrained to lie beside its root node.

        :param iel: ideal edge length for the graph
        :return: nothing
        """
        self.face.insertTreeNode(self, iel)


class Side:

    def __init__(self, face, nodeseq, direc):
        # the Face to which this Side belongs:
        self.face = face
        # the Nodes belonging to this Side, in clockwise order w.r.t. the Face:
        self.nodeseq = nodeseq
        # the Compass direc in which we move as we go forward through the nodeseq:
        self.forward = direc
        # the Compass direc pointing into the interior of the Face:
        self.inward = Compass.cw90(direc)
        # the varying dimension:
        self.vardim = Compass.variableDimension[direc]
        # the constant dimension:
        self.constdim = Compass.constantDimension[direc]
        # a list in which to store TreePlacement objects when trees are placed on
        # this Side:
        self.treePlacements = []

    def getNumRootNodes(self):
        return len(filter(lambda v: v.isRootNode, self.nodeseq))

    def containsNode(self, node):
        """
        :param node: a Node object
        :return: boolean saying if that Node belongs to this Side
        """
        return node in self.nodeseq

    def containsSubseq(self, subseq):
        """
        :param subseq: list of nodes
        :return: boolean, saying whether the nodeseq of this Side contains the given subseq as
                 an ordered subsequence
        """
        if len(subseq) == 0: return True
        u0 = subseq[0]
        if not u0 in self.nodeseq: return False
        i0 = self.nodeseq.index(u0)
        n = len(subseq)
        S = self.nodeseq[i0:i0+n]
        return S == subseq

    def firstNode(self):
        return self.nodeseq[0]

    def aWidestNode(self):
        """
        :return: any Node on this side which is widest among all the nodes on this side,
                 "width" being in the dimension transverse to the orientation of the side
        """
        width = (lambda n: n.w) if self.forward in Compass.vertical else (lambda n: n.h)
        return max(self.nodeseq, key=width)

    def constCoord(self):
        """
        :return: the constant coord for this side, shared by all edges

        NB: For now we are working under the assumption that all edges connect
        at the /centres/ of the sides of all nodes!
        """
        n0 = self.nodeseq[0]
        z = n0.x if self.constdim == XDIM else n0.y
        return z

    def liesOppositeSegment(self, seg, openInterval=False):
        """
        :param seg: a LineSegment
        :return: boolean saying whether the closed interval spanned by this Side
                 runs in the same dimension as the segment, and overlaps it in
                 projection onto that dimension
        """
        return self.getIntervalOppositeSegment(seg, openInterval=openInterval) is not None

    def getFirstPtOppositeSegment(self, seg):
        """
        :param seg: a LineSegment
        :return: the first point of the interval of this Side that lies opposite
                 the given segment, if any, else None
        """
        I = self.getIntervalOppositeSegment(seg)
        if I is None:
            pt = None
        else:
            w = I[0]
            z = self.constCoord()
            pt = (z, w) if self.constdim == XDIM else (w, z)
        return pt

    def getIntervalOppositeSegment(self, seg, openInterval=False):
        """
        :param seg: a LineSegment
        :return: interval I = [a, b] being the intersection of this Side's
                 closed interval with that of the given segment, or None if the
                 intersection is empty or if this Side is not even parallel to
                 the segment
        """
        # First establish whether we run parallel.
        if not Compass.sameDimension(self.forward, seg.direc): return None
        # Now consider where the intervals overlap.
        if openInterval:
            I = seg.openIntervalIntersection(self.closedInterval())
        else:
            I = seg.closedIntervalIntersection(self.closedInterval())
        return I

    def closedInterval(self):
        """
        :return: a closed interval I = [a, b] where a and b are the extreme
        coordinates covered by this Side, up to the extremes of the boxes of the
        extreme nodes
        """
        x, X, y, Y = self.nodeseq[0].boundingBoxxXyY()
        u, U, v, V = self.nodeseq[-1].boundingBoxxXyY()
        if self.forward == Compass.EAST:
            I = [x, U]
        elif self.forward == Compass.SOUTH:
            I = [y, V]
        elif self.forward == Compass.WEST:
            I = [u, X]
        else:
            assert(self.forward == Compass.NORTH)
            I = [v, Y]
        return I

    def halfwidthOppositeSegment(self, sign, seg):
        """
        :param sign: an integer +/-1 indicating whether we are interested in the
                     upper or lower half-width, respectively
        :param seg: a LineSegment
        :return: the half-width of this Side where it overlaps the given segment,
                 or None if there is no overlap

        Note: for now, while edges always meet the centres of node sides, the sign
        argument is irrelevant.
        """
        I = self.getIntervalOppositeSegment(seg)
        # If we are not opposite the segment, return None.
        if I is None: return None
        a, b = I
        # Set half the thickness of an edgenode as the minimum halfwidth.
        hw = Graph.EDGENODE_THICKNESS / 2.0
        # Prepare functions for reading data off of nodes.
        if self.forward in Compass.horizontal:
            interval = (lambda t: t[:2])
            dimension = (lambda u: u.w)
        else:
            interval = (lambda t: t[2:])
            dimension = (lambda u: u.h)
        # Now consider all nodes belonging to the Side.
        for node in self.nodeseq:
            c, d = interval(node.boundingBoxxXyY())
            if b >= c and d >= a:
                # If the node is present in the interval I, consider its size.
                hw = max(hw, dimension(node)/2.0)
        return hw

    def __iter__(self):
        return iter(self.nodeseq)

    def treePlacementIsRelevant(self, tp):
        """
        :param tp: a TreePlacement whose root node belongs to this Side
        :return: boolean saying if the placement is actually relevant to this
                 side

                 This will always be True if the root node is not a pinch
                 node, and there are no bridge edges, but otherwise it may be False.
        """
        root = tp.node
        assert(root in self.nodeseq)
        pairs = self.face.nbrPairs[root.ID]
        if len(pairs) == 1:
            return True
        else:
            assert(len(pairs) == 2)
            # In this case there are four Sides of self.Face which include root
            # among their nodes, and self is among those.
            # But the tree placement tp is relevant to just two of those Sides,
            # and we must decide if self is one of them.
            # We consider the four nodes named in the neighbour pairs.
            a, b = pairs[0]
            c, d = pairs[1]
            nbrs = [a, b, c, d]
            # The order of the set {a, b, c, d} may be four, but it may also be
            # two or three, if root has one resp. two bridge edges.
            # Nevertheless, since the neighbour pairs were computed by traversing the
            # nodes of the face in the forward order, the list [a, b, c, d] is in
            # anticlockwise order, where we understand the outgoing side of a bridge
            # edge to come before the incoming side, in anticlockwise order.
            # With this understanding, we simply need to compute the compass direction
            # from root to each of a, b, c, d, and then find where the placement direction
            # fits into this order. The two adjacent sides are the two to which the tree placement
            # is relevant.
            direcs = [Compass.cardinalDirection(root, v) for v in nbrs]
            dp = tp.placementDirec
            # We should not be considering a tree placement in a direction already occupied
            # by one of the sides of the face.
            assert(dp not in direcs)
            # Now compute the indices i0, i1 in {0, 1, 2, 3} of the two nbrs which are the closest
            # to dp in each direction.
            # For a given direction we cannot use the mere compass rose distance, since that does
            # not take account of whether the edge is outgoing or incoming.
            # So we double the rose distance and add 1 for the "more distant direction," which
            # is the outgoing edge when looking clockwise, and the incoming edge when looking in
            # the anticlockwise direction.
            # We can judge incoming vs. outgoing by the parity of index i, since nbr pairs are
            # always of the form [incoming, outgoing].
            i0 = min(range(4),
                key=lambda i: 2*Compass.cwRoseDistance(dp, direcs[i]) + (1 if i%2==1 else 0)
            )
            i1 = min(range(4),
                key=lambda i: 2*Compass.acwRoseDistance(dp, direcs[i]) + (1 if i%2==0 else 0)
            )
            # Finally we must say whether self is among the two Sides corresponding to i0, i1.
            for i in [i0, i1]:
                subseq = [nbrs[i], root]
                if i%2 == 1: subseq.reverse()
                if self.containsSubseq(subseq):
                    # In this case self is one of the relevant sides.
                    return True
            else:
                # self was not one of the relevant sides.
                return False

    def addTreePlacement(self, tp):
        assert(tp.node in self.nodeseq)
        self.treePlacements.append(tp)

    def computeCollateralProjSeq(self, tp0, iel):
        """
        :param tp0: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: a ProjSeq to remove/prevent overlaps between the given TreePlacement's
                 tree box, and any existing treenodes already on this Side, as well as
                 ordinary perimeter nodes on the Side.

                 Replaces old method, computeCollateralTreeSep.
        """
        assert(tp0.node in self.nodeseq)
        ps = ProjSeq()
        tw, th, tu, tv = tp0.treeBoxWithRootVector(iel=iel)
        root0 = tp0.node
        rw, rh, rx, ry = root0.w, root0.h, root0.x, root0.y
        # Compute segment representing extent of new treenode in Side's inward direction.
        if self.inward == Compass.NORTH:
            seg = LineSegment((rx, ry-rh/2.0), (rx, ry+tv-th/2.0))
        elif self.inward == Compass.EAST:
            seg = LineSegment((rx+rw/2.0, ry), (rx+tu+tw/2.0, ry))
        elif self.inward == Compass.SOUTH:
            seg = LineSegment((rx, ry+rh/2.0), (rx, ry+tv+th/2.0))
        elif self.inward == Compass.WEST:
            seg = LineSegment((rx-rw/2.0, ry), (rx+tu-tw/2.0, ry))
        # May need constraints only if the segment points inward w.r.t. this Side.
        if seg.direc == self.inward:
            pcs = []
            # Handle nodes on each side of root node.
            i0 = self.nodeseq.index(root0)
            # Perimeter nodes:
            before, after = self.nodeseq[:i0], self.nodeseq[i0+1:]
            # Treenodes:
            for tp1 in self.treePlacements:
                i1 = self.nodeseq.index(tp1.node)
                assert(tp1.boxNode is not None)
                assert(i1 != i0)
                if i1 < i0:
                    before.append(tp1.boxNode)
                elif i1 > i0:
                    after.append(tp1.boxNode)
            # Filter out nodes that do not lie opposite the segment computed above.
            before = filter(lambda n: n.liesOppositeSegment(seg, openInterval=True), before)
            after = filter(lambda n: n.liesOppositeSegment(seg, openInterval=True), after)
            # Prepare "signs" to manage the work in a for-loop.
            signs = (-1, 1) if self.forward in Compass.increasing else (1, -1)
            jobs = zip(signs, (before, after))
            for sign, nodes in jobs:
                for node in nodes:
                    left, right = (root0, node) if sign == 1 else (node, root0)
                    # Compute gap.
                    if self.forward in Compass.horizontal:
                        gap = node.w/2.0 + sign*tu + tw/2.0
                    else:
                        gap = node.h/2.0 + sign*tv + th/2.0
                    pcs.append(SepCo(self.vardim, left, right, gap))
            if len(pcs) > 0:
                ps.addConstraintSet(pcs, self.vardim)
        return ps

class BoundaryCrossing:

    def __init__(self, x, y, owners):
        self.x = x
        self.y = y
        self.owners = owners

class BoundarySegment:

    def __init__(self, p0, cross0, p1, cross1, owners, degenerateDirec=None):
        """
        :param p0: coords (x0, y0) of first endpt of the segment
        :param cross0: boolean saying if we cross the face boundary at p0
        :param p1: coords (x1, y1) of second endpt of the segment; should be
                   clockwise from p0, in traversing the boundary
        :param cross1: boolean saying if we cross the face boundary at p1
        :param owners: list of the one or two nodes associated with this segment:
                       one if the segment is part of a node boundary, two if it
                       represents an edge
        :param degenerateDirec: Set this to a cardinal Compass direction if the
         segment is degenerate in that its two endpoints are the same. This is
         provided for cases where we have such degenerate segments and still want
         to think of them as having a direction.
        """
        # Either we need distinct points p0, p1, or else you have to have
        # said what the direction is to be, for a degenerate segment.
        assert(p0 != p1 or degenerateDirec is not None)
        EPSILON = 0.00001
        self.owners = owners
        x0, y0 = p0
        x1, y1 = p1
        if p0 == p1:
            dd = degenerateDirec
            self.direc = dd
            lowCross, highCross = cross0, cross0
            if dd in Compass.vertical:
                self.constCoord = x0
                self.lowEnd, self.highEnd = y0, y0
            else:
                self.constCoord = y0
                self.lowEnd, self.highEnd = x0, x0
        else:
            if abs(x0 - x1) < EPSILON:
                self.constCoord = x0
                if y0 < y1:
                    self.direc = Compass.SOUTH
                    self.lowEnd, self.highEnd = y0, y1
                    lowCross, highCross = cross0, cross1
                else:
                    assert(y1 < y0)
                    self.direc = Compass.NORTH
                    self.lowEnd, self.highEnd = y1, y0
                    lowCross, highCross = cross1, cross0
            else:
                assert(abs(y0 - y1) < EPSILON)
                self.constCoord = y0
                if x0 < x1:
                    self.direc = Compass.EAST
                    self.lowEnd, self.highEnd = x0, x1
                    lowCross, highCross = cross0, cross1
                else:
                    assert(x1 < x0)
                    self.direc = Compass.WEST
                    self.lowEnd, self.highEnd = x1, x0
                    lowCross, highCross = cross1, cross0
        self.crossings = {-1: lowCross, 0: True, 1: highCross}

    def __repr__(self):
        s = ''
        s += 'Seg: %s from %.1f to %.1f at %.1f, with crossings: %s, %s' % (
            '|' if self.direc in Compass.vertical else '--',
            self.lowEnd, self.highEnd, self.constCoord,
            self.crossings[-1], self.crossings[1]
        )
        return s

    def __lt__(self, other):
        return self.constCoord < other.constCoord

    def isTransverseTo(self, direc):
        """
        :param direc: a cardinal Compass direction
        :return: boolean saying if this segment is transverse to that direction, i.e.
                 if this segment's direction is in the complementary dimension
        """
        return not Compass.sameDimension(direc, self.direc)

    def closedIntervalContains(self, w):
        """
        :param w: a float
        :return: say whether the closed interval contains w
        """
        return self.lowEnd <= w and w <= self.highEnd

    def openIntervalContains(self, w):
        return self.lowEnd < w and w < self.highEnd

    def crossesOpenInterval(self, z0, z1):
        """
        :param z0: float
        :param z1: float
        :return: boolean saying if the constant coord lies in the open interval
        """
        return z0 < self.constCoord and self.constCoord < z1

    def getCrossingAt(self, w):
        assert(self.closedIntervalContains(w))
        z = self.constCoord
        if self.direc in Compass.horizontal:
            return (w, z)
        else:
            return (z, w)

    def getConstCoord(self):
        return self.constCoord

    def canCrossBoundaryAt(self, w):
        """
        :param w: float
        :return: boolean saying whether the boundary of the face is crossed at w
        """
        if not self.closedIntervalContains(w): return False
        if w == self.lowEnd:
            i = -1
        elif w == self.highEnd:
            i = 1
        else:
            assert(self.openIntervalContains(w))
            i = 0
        return self.crossings[i]


class Face:

    def __init__(self, G, logger, config):
        self.graph = G
        self.logger = logger
        self.config = config
        self.nodeseq = []
        self.n = 0
        self.external = False
        self.ID = None

        # Dictionary nodeID:list
        # list is of the form [(w, v)] or [(w1, v1), (w2, v2)]
        # giving the one or two (for pinch nodes) neighbour pairs of the
        # key node, in clockwise order.
        self.nbrPairs = {}

        self.sides = []
        self.nodeIDsToSides = {}

        # Dictionary of Nodes, by ID, being the large boxes added to represent
        # whole trees, with padding:
        self.treenodes = {}

        self.nodeIDToTreePlacement = {}

    def getAllTreePlacements(self):
        return self.nodeIDToTreePlacement.values()

    def applyGeometryToTrees(self, iel=0):
        for tp in self.nodeIDToTreePlacement.values():
            tp.applyGeometryToTree(iel=iel)

    def insertTreesIntoGraph(self, H):
        """
        :param H: a Graph
        :return: nothing

        Add the individual nodes of the trees to the graph, and remove the
        corresponding treenode.
        """
        for tp in self.nodeIDToTreePlacement.values():
            tp.insertTreeIntoGraph(H)

    def numTreesWithGrowthDirec(self, dg, scaleBySize=False):
        tps = self.nodeIDToTreePlacement.values()
        tps = filter(lambda tp: tp.growthDirec == dg, tps)
        if scaleBySize:
            return sum(tp.tree.size() for tp in tps)
        else:
            return len(tps)

    def getAllTreePlacements(self):
        return self.nodeIDToTreePlacement.values()

    def getAllSidesOppositeSegment(self, seg, openInterval=False):
        return filter(lambda s: s.liesOppositeSegment(seg, openInterval=openInterval), self.sides)

    def getAllTreenodesOppositeSegment(self, seg, openInterval=False):
        return filter(lambda t: t.liesOppositeSegment(seg, openInterval=openInterval), self.treenodes.values())

    def getAllPerimeterNodesOppositeSegment(self, seg, openInterval=False):
        return filter(lambda n: n.liesOppositeSegment(seg, openInterval=openInterval), self.nodeseq)

    def buildOrthoRoutingRigForSolidifiedFace(self, iel):
        # To ensure that there are alleys between treenodes and edgenodes,
        # here we make the edgenodes /half/ as thick as they are when we
        # are shaking the graph with solidified edges.
        EDGENODE_THICKNESS = Graph.EDGENODE_THICKNESS / 2.0
        rr = RoutingRig({
            'routing': OrthogonalRouting
        })
        nodes = {}
        maxID = -1
        for node in self.nodeseq:
            nodes[node.ID] = node
            maxID = max(maxID, node.ID)
        nextID = maxID + 1
        for node in self.nodeseq:
            if node.ID in self.nodeIDToTreePlacement:
                tp = self.nodeIDToTreePlacement[node.ID]
                w, h, u, v = tp.treeBoxWithRootVector(iel=0)
                dummy = Node()
                dummy.x, dummy.y = node.x + u, node.y + v
                dummy.w, dummy.h = w, h
                dummy.ID = nextID
                nextID += 1
                nodes[dummy.ID] = dummy
        maxID = nextID - 1

        baseID = maxID + 1
        for i in range(self.n):
            u, v = self.nodeseq[i], self.nodeseq[(i + 1) % self.n]
            node = Node()
            node.ID = baseID + i
            direc = Compass.cardinalDirection(u, v)
            # Unlike in the Graph.solidifyEdges method, we will make these edgenodes
            # overlap their endpoint nodes. This is necessary so that the router
            # will not try to thread a path in between.
            # These nodes are only going to be used to build ShapeRefs for the
            # router, so their overlapping just leaves the router with one giant
            # compound shape representing the entire face with nodes and edges.
            if direc in Compass.vertical:
                node.w = EDGENODE_THICKNESS
                node.h = abs(u.y - v.y)
            else:
                node.h = EDGENODE_THICKNESS
                node.w = abs(u.x - v.x)
            node.x = (u.x + v.x) / 2.0
            node.y = (u.y + v.y) / 2.0
            nodes[node.ID] = node
        rr.addNodesAndEdges(nodes, {})
        return rr

    def initRoutingRig(self, routerOpts, fixedStraightRoutes=False):
        """
        Set up a RoutingRig object for the nodes and edges of this face.
        :param routerOpts: dictionary of options for the router
        :param fixedStraightRoutes: set True if you want the edges to have fixed
                                    straight routes
        :return: a RoutingRig object
        """
        rr = RoutingRig(routerOpts)
        nodes = {}
        for node in self.nodeseq:
            nodes[node.ID] = node
        edges = {}
        for i in range(self.n):
            u, v = self.nodeseq[i], self.nodeseq[(i + 1) % self.n]
            edge = self.graph.getEdgeBtwNodes(u, v)
            assert(edge is not None)
            edges[repr(edge)] = edge
        rr.addNodesAndEdges(nodes, edges, fixedStraightRoutes=fixedStraightRoutes)
        return rr

    def direc(self, u, v):
        return self.graph.nodeConf.getDirec(u, v)

    def initWithEdgeSeq(self, edges):
        """
        :param edges: a list of ordered pairs
                 (s0, t0), (s1, t1), ..., (sn-1, tn-1)
            of node IDs, such that tk = sk+1 for all k from 0 to n-1 inclusive,
            where k+1 is understood mod n
        :return: nothing
        """
        self.nodeseq = [self.graph.nodes[e[0]] for e in edges]
        self.n = len(self.nodeseq)
        self.computeNbrPairs()
        self.computeSides()

    def computeNbrPairs(self):
        # First pass: identify the one or two indices at which each node
        # occurs in the nodeseq.
        indices = {}
        for i, v in enumerate(self.nodeseq):
            L = indices.get(v.ID, [])
            L.append(i)
            indices[v.ID] = L
        # Now assemble the neighbour pairs.
        n = self.n
        ns = self.nodeseq
        for v in ns:
            L = indices[v.ID]
            pairs = [(ns[(i-1)%n], ns[(i+1)%n]) for i in L]
            self.nbrPairs[v.ID] = pairs

    def computeSides(self):
        # Get index of first bend node.
        for i0 in range(self.n):
            if self.isBendNode(self.nodeseq[i0]):
                break
        # Prepare a node sequence starting from node i0, then
        # cycling back around to i0 and including node i0 /again/.
        # Then a Side begins at the start of this sequence, and a
        # Side ends at the end of it, and there is no Side between
        # the last and first one.
        ns = self.nodeseq[i0:] + self.nodeseq[:i0+1]
        # Set up the loop.
        assert(self.n >= 3)
        u, v = ns[:2]
        nodes = [u, v]
        d0 = self.direc(u, v)
        ns = ns[2:]
        for v in ns:
            u = nodes[-1]
            d1 = self.direc(u, v)
            if d1 == d0:
                nodes.append(v)
            else:
                self.sides.append(Side(self, nodes, d0))
                nodes = [u, v]
                d0 = d1
        # Create the final Side.
        self.sides.append(Side(self, nodes, d0))
        # Now index the Sides by the Nodes they contain.
        for S in self.sides:
            for u in S:
                L = self.nodeIDsToSides.get(u.ID, [])
                L.append(S)
                self.nodeIDsToSides[u.ID] = L

    def __repr__(self):
        return 'Face%s: %s' % (
            ' (external)' if self.external else '',
            [node.ID for node in self.nodeseq]
        )

    def isBendNode(self, v):
        pairs = self.nbrPairs[v.ID]
        if len(pairs) == 2:
            return True
        else:
            u, w = pairs[0]
            duv = self.graph.nodeConf.getDirec(u, v)
            dvw = self.graph.nodeConf.getDirec(v, w)
            return duv != dvw

    def computeBoundarySegments(self, iel):
        segs = []
        # Since a node may occur twice in the traversal of a face, we keep a list
        # of nodes visited, to prevent duplicate boundary segments (which wreak all
        # sorts of havoc -- e.g. inside/outside is judged by parity of number of
        # crossings, and this goes haywire if a crossing is counted twice!).
        visitedIDs = set()
        for u in self.nodeseq:
            if u.ID in visitedIDs: continue
            visitedIDs.add(u.ID)
            pairs = self.nbrPairs[u.ID]
            for w, v in pairs:

                # Begin with the edge from u to v.
                outdir = self.direc(u, v)
                p0, p1 = u.getBdryCompassPt(outdir), v.getBdryCompassPt(Compass.flip(outdir))
                # If nodes u and v are touching, then we will get a degenerate segment.
                # This is still usable, but we have to tell it what its "direction" is.
                segs.append(BoundarySegment(p0, True, p1, True, [u, v], degenerateDirec=outdir))

                # Now do the segments for node u traversed when proceeding anticlockwise from
                # the edge connecting to neighbour w, to the edge connecting to neighbour v.
                #
                # WARNING: this code depends on the assumption that edges always
                # attach to nodes at midpoints of the sides of their bounding boxes.
                #
                if u.ID in self.nodeIDToTreePlacement:
                    tp = self.nodeIDToTreePlacement[u.ID]
                    # We use the /unpadded/ treenodes when considering where the boundary lies:
                    padding = 0
                    segs.extend(tp.computeBoundarySegments(padding))
                else:
                    # Get the direction from u back to the previous node.
                    backdir = self.direc(u, w)
                    assert(backdir != outdir)
                    # Get the set of all compass directions strictly between (non-inclusive)
                    # the direction to w and the direction to v, going anticlockwise.
                    # NB: We go anticlockwise around the node so that the /right-hand/
                    # side of its boundary segments will be the interior of the face.
                    rr = Compass.acwRose + Compass.acwRose
                    i0 = rr.index(backdir)
                    i1 = i0 + rr[i0:].index(outdir)
                    btw = rr[i0+1:i1]
                    # Now the directions we want are that to w, then all and only the /ordinal/
                    # directions between, and then the direction to v.
                    # Throw any cardinal directions out of the btw list.
                    ordsbtw = filter(lambda d: d in Compass.cwOrds, btw)
                    # And now tack the directions to v and to w onto the ends of the list.
                    dirs = [backdir] + ordsbtw + [outdir]
                    # The interior segments of the node are those running between the
                    # boundary points in consecutive directions in this list of dirs.
                    pts = [u.getBdryCompassPt(d) for d in dirs]
                    for j in range(len(pts) - 1):
                        p0, p1 = pts[j:j+2]
                        d0, d1 = dirs[j:j+2]
                        # We consider node boundary segments to be interior, but edge segments
                        # to be exterior (because things work well this way), which means
                        # there is no crossing at either end of a node boundary segment.
                        segs.append(BoundarySegment(p0, False, p1, False, [u]))
        return segs


    def insertTreeNode(self, tp, iel):
        """
        To be used after the face has been expanded to make room for the tree.
        This method adds a large node to the graph, representing the bounding
        box of the tree, with padding.
        The treenode is constrained to lie beside its root node.

        :param tp: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: nothing
        """
        # Ask Side object(s) to note the placement.
        for S in self.nodeIDsToSides[tp.node.ID]:
            if S.treePlacementIsRelevant(tp):
                S.addTreePlacement(tp)
        # Now insert the treenode.
        w, h, u, v = tp.treeBoxWithRootVector(iel=iel)
        # Create the treenode and add it to the graph.
        root = tp.node
        treenode = Node()
        treenode.w, treenode.h = w, h
        treenode.x, treenode.y = root.x + u, root.y + v
        treenode.ID = self.graph.getNextID()
        treenode.fill = '#C0804080'  # alpha = 0.5 for transparency
        treenode.setIDAsLabel()
        self.graph.addNode(treenode)
        # Make records
        self.treenodes[treenode.ID] = treenode
        tp.boxNode = treenode
        self.nodeIDToTreePlacement[root.ID] = tp
        # Constrain the treenode to sit beside the root node.
        d = tp.growthDirec
        offset = u if d in Compass.vertical else v
        self.graph.nodeConf.setDirec(root, treenode, d + 8,
                                     alignOffsets=(0, offset))

    def applyPS(self, ps, iel):
        """
        :param ps: a ProjSeq object
        :param iel: ideal edge length
        :return: boolean saying whether all the projections were successful

        Convenience function for applying a projseq with all the necessary options
        """
        return self.graph.applyProjSeq(ps, self.logger, iel=iel, opx=True, opy=True,
                                    solidEdgesX=True, solidEdgesY=True,
                                    accept=0)

    def evaluateExpansionOptions(self, tp, iel):
        """
        :param tp: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: list of ProjSeq objects, with stress costs
        """
        logger = self.logger
        # Save initial node positions.
        self.graph.pushNodePoses()
        # Start by removing any overlaps with collateral treenodes.
        ps0 = self.doCollateralExpansion(tp, iel)
        # Now consider the options for removing remaining overlaps.
        projseqs = []
        if self.config.HEURISTIC_CHOICE_FOR_PRIMARY_EXPANSION_DIMENSION:
            wd = WaterDivide(tp, XDIM, iel, logger)
            xEst, yEst = wd.estimateCostByDimension()
            if self.config.HCPED_COSTLIER_DIMENSION_FIRST:
                # We work in the coslier dimension first, hoping that the bigger
                # change that it makes will already handle the other dimension.
                dims = [XDIM] if xEst > yEst else [YDIM]
            else:
                # We work in the cheaper dimension first, hoping it's enough and
                # we get away with it.
                dims = [XDIM] if xEst < yEst else [YDIM]
        else:
            # We'll try each dimension as the initial dimension.
            dims = [XDIM, YDIM]
        self.graph.pushNodePoses()
        if logger.level >= LogLevel.TIMING:
            suffix = logger.lastSuffix
        for dim in dims:
            self.graph.popNodePoses()
            self.graph.pushNodePoses()
            if logger.level >= LogLevel.TIMING:
                logger.lastSuffix = suffix
            wd = WaterDivide(tp, dim, iel, logger)
            projseqs.extend(wd.buildProjSeqs(ps0))
        self.graph.dropNodePoses()
        # Restore node positions and return.
        self.graph.popNodePoses()
        return projseqs

    def estimateExpansionOptions(self, tp, iel):
        """
        Like the evaluateExpansionOptions method, but only estimates costs,
        instead of carrying out all projections and backtracking.
        :param tp: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: list of ProjSeq objects, with stress costs
        """
        logger = self.logger
        # Start by removing any overlaps with collateral treenodes.
        ps0 = self.doCollateralExpansion(tp, iel, testFeasibility=False)
        # Now consider the options for removing remaining overlaps.
        projseqs = []
        if self.config.HEURISTIC_CHOICE_FOR_PRIMARY_EXPANSION_DIMENSION:
            wd = WaterDivide(tp, XDIM, iel, logger)
            xEst, yEst = wd.estimateCostByDimension()
            if self.config.HCPED_COSTLIER_DIMENSION_FIRST:
                # We work in the coslier dimension first, hoping that the bigger
                # change that it makes will already handle the other dimension.
                dims = [XDIM] if xEst > yEst else [YDIM]
            else:
                # We work in the cheaper dimension first, hoping it's enough and
                # we get away with it.
                dims = [XDIM] if xEst < yEst else [YDIM]
        else:
            # We'll try each dimension as the initial dimension.
            dims = [XDIM, YDIM]
        if logger.level >= LogLevel.TIMING:
            suffix = logger.lastSuffix
        for dim in dims:
            if logger.level >= LogLevel.TIMING:
                logger.lastSuffix = suffix
            wd = WaterDivide(tp, dim, iel, logger)
            projseqs.extend(wd.buildProjSeqs(ps0, estimate=True))
        # Now we must discard any infeasible projseqs at the beginning of the list;
        # i.e. at least the first one in the list must be feasible.
        def firstIsFeasible():
            self.graph.pushNodePoses()
            okay = self.graph.applyProjSeq(projseqs[0], logger, opx=True, opy=True, solidEdgesX=True, solidEdgesY=True)
            self.graph.popNodePoses()
        while not firstIsFeasible():
            projseqs = projseqs[1:]
        return projseqs

    def doCollateralExpansion(self, tp, iel, testFeasibility=True):
        """
        :param tp: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: an optimal ProjSeq with which to remove/prevent overlap
                 with nodes and treenodes belonging to the same side or sides
                 to which the root node of the TreePlacement belongs.

                 The returned ProjSeq has already been evaluated and applied.

                 Raises an exception if any of the ProjSeqs attempted is
                 infeasible.
        """
        logger = self.logger

        def tryProjSeq(ps, idx):
            if logger.level >= LogLevel.TIMING:
                suffix = '_exp_N%s_P%s_G%s_%s_col%d' % (
                    tp.node.ID, tp.placementDirec, tp.growthDirec,
                    'fl' if tp.flip else 'nf',
                    idx
                )
                if logger.level >= LogLevel.PROGRESS:
                    print suffix
                logger.startNewTimer('projection')
            ok = self.applyPS(ps, iel)
            if logger.level >= LogLevel.TIMING:
                logger.stopLastTimer()
            if logger.level >= LogLevel.FINER_STAGE_GRAPHS:
                logger.writeGML(suffix)
            if not ok:
                raise Exception('Infeasible collateral tree sep!')

        PS = self.getCollateralProjSeqs(tp, iel)
        if len(PS) == 1:
            ps0 = PS[0]
            if testFeasibility:
                tryProjSeq(ps0, 0)
        else:
            raise Exception('Did not get unique collateral projseq.')
        return ps0

    def getRelevantSidesForPlacement(self, tp):
        """
        :param tp: a TreePlacement
        :return: a list of all the Sides that are relevant for this TreePlacement
        """
        sides = self.nodeIDsToSides[tp.node.ID]
        sides = filter(lambda S: S.treePlacementIsRelevant(tp), sides)
        return sides

    def getCollateralProjSeqs(self, tp, iel):
        """
        :param tp: a TreePlacement object
        :param iel: ideal edge length for the graph
        :return: a list of ProjSeqs to remove overlaps with pre-existing treenodes

        NB: In the case where tp's root lies on two sides of the face, these two
        sides must be aligned in complementary dimensions (i.e. one in x, and one
        in y), so the sepcos generated here are always independent.

        That is, if each side has an existing treenode that must be pushed away,
        the sepco for one tree will never achieve the push for the other, since
        these pushes /must/ be done in complementary dimensions.

        Therefore there is no need to wait until after the first removal has been
        actually performed to compute the sepco for the second removal; it cannot
        be affected by the first one.

        I don't think it matters the order in which we do these two projections either.
        """
        pses = []
        sides = self.getRelevantSidesForPlacement(tp)
        assert(len(sides) in [1, 2])
        if len(sides) == 1:
            pses.append(sides[0].computeCollateralProjSeq(tp, iel))
        else:
            Q = [S.computeCollateralProjSeq(tp, iel) for S in sides]
            ps = ProjSeq()
            for qs in Q:
                ps += qs
            pses.append(ps)
        return pses

    def inwardDirsAvail(self, node):
        """
        :param node: a Node belonging to the face
        :return: a list of Compass directions in which an edge could point
                 if it were anchored at the given node and pointed inward,
                 into the face

                 Note: if any cardinal direction(s) are available then we
                 will NOT report available ordinal directions (diagonal)
        """
        dirs = []
        pairs = self.nbrPairs[node.ID]
        for pair in pairs:
            prev, next = pair
            # Because of the way the faces are traced out, swinging from the previous
            # node to the next one sweeps out an area of the /inside/ of this face
            # when you go /anticlockwise/.
            #
            # Therefore if D0 is the compass direction from node to prev, and if D1
            # is that from node to next, then the available directions are all those
            # strictly between D0 and D1 /in the anticlockwise direction/.
            #
            # For example, if node --> prev points South while node --> next points
            # West, then both East and North are available.
            # But if node --> prev points East while node --> next points North, then
            # only NE is available.
            D0 = Compass.cardinalDirection(node, prev)
            D1 = Compass.cardinalDirection(node, next)
            rr = Compass.acwRose + Compass.acwRose
            i0 = rr.index(D0)
            i1 = i0 + rr[i0:].index(D1)
            btw = rr[i0+1:i1]
            if set(btw).isdisjoint(set(Compass.cwCards)):
                # btw contains NO cardinal directions
                pass
            else:
                # btw contains at least one cardinal direction,
                # so take away the ordinal directions.
                btw = list(set(btw) - set(Compass.cwOrds))
            dirs.extend(btw)
        return dirs

    def allPossibleTreePlacements(self, tree, node):
        """
        :param tree: a Tree to be placed
        :param node: a Node belonging to the face
        :return: a list of TreePlacement objects giving all the possible placements
                 of a tree into this face at the given node, but with cost not yet
                 evaluated
        """
        dirs = self.inwardDirsAvail(node)
        tps = []
        for d1 in dirs:
            if d1 in Compass.cwCards:
                tps.append(TreePlacement(tree, self, node, d1, d1))
                if not tree.isSymmetric:
                    tps.append(TreePlacement(tree, self, node, d1, d1, flip=True))
            else:
                for d2 in Compass.components[d1]:
                    tps.append(TreePlacement(tree, self, node, d1, d2))
                    if not tree.isSymmetric:
                        tps.append(TreePlacement(tree, self, node, d1, d2, flip=True))
        return tps
