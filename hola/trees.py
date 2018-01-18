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
from collections import deque as Deque
from ortho import *
from constraints import AlignCo, SepCo, DistCo, FlexDistCo, FixedRelCo
from graphs import Node
from edgeoverlaps import tolerantPartition

class Tree:

    def __init__(self, G, r):
        # G must be an acyclic graph or we will raise an exception.
        # (At least it must be acyclic south of the nominated root.)
        self.graph = G
        self.root = r
        self.nodes = {} # ID of node --> node
        # Initialise fields.
        self.depth = 0
        self.breadth = 0
        self.leaves = {} # ID of node --> 1
        self.nodesByRank = {} # rank (int) --> list of Nodes
        self.rankByNodeID = {} # ID of node --> rank (int)
        self.parent = {} # ID of child --> parent Node
        self.isSymmetric = False
        self.holaConfig = None
        # For layout, we have a map
        # rank (int) --> [lb, ub]
        # returning the lower and upper bounds on the lateral coordinates
        # of the tree, for each rank (e.g. for NORTH growth direction the
        # bounds are on x-coordinates).
        self.boundsByRank = {}
        # We also keep the global lower and upper bounds over all populated ranks.
        self.lb = 0
        self.ub = 0
        self.orderingNumbers = {} # ID of node --> its ordering number
        self.bufferNodes = []
        self.pcs = []
        # The following booleans configure how the getBounds function for this
        # tree will work.
        #
        # If the boundary is infinite, then you will always get
        # lower and upper bound for every rank, no matter how high; else you will
        # get None for ranks having no nodes in them.
        #
        # If the boundary is tight then the bounds for each rank cover just the
        # nodes on that rank plus half nodeSep on each end; else the bounds for
        # every rank are equal to the tight bounds for the widest rank.
        #
        # The style of the layout can be configured to some extent by setting
        # these booleans.
        self.boundaryInfinite = False
        self.boundaryTight = True
        #
        self.growthDir = None
        # Compute ranks etc.
        self.setRank(r,0)
        #queue = [r]
        deque = Deque([r])
        #
        while deque:
            # Pop front of queue.
            node = deque.popleft()
            # Make sure we're not cycling.
            if node.ID in self.nodes:
                msg = 'Trying to construct a Tree on a '
                msg += 'graph with a cycle.'
                raise Exception(msg)
            self.nodes[node.ID] = node
            # Proceed.            
            children = node.getChildren()
            rank = self.rankByNodeID[node.ID]
            if not children:
                self.leaves[node.ID] = node
            for c in children:
                self.parent[c.ID] = node
                self.setRank(c, rank+1)
                deque.append(c)
        ranks = self.nodesByRank.keys()
        self.depth = len(ranks)
        self.breadth = max([
            len(r) for r in self.nodesByRank.values()
        ])

    def setHolaConfig(self, hc):
        self.holaConfig = hc

    def cmpRankmates(self, u, v):
        """
        If u and v are two nodes belonging to the same rank in this tree,
        then return -1, 0, 1 to indicate their order, so that the rank can
        be sorted. This ordering is based on the ordering numbers that should
        have been assigned by now, by the computeOrderingNumbers function.
        """
        pu, pv = self.parent[u.ID], self.parent[v.ID]
        if pu == pv:
            return self.orderingNumbers[u.ID] - self.orderingNumbers[v.ID]
        else:
            return self.cmpRankmates(pu, pv)

    def computeAlignedSets(self, root):
        """
        :param root: a Node
        :return: dictionary {'h': Hsets, 'v': Vsets}, where Hsets is a list of lists,
                 each memeber of which is a list of nodes that are horizontally aligned;
                 likewise with Vsets.

                 The supplied root is substituted for this tree's root node.
        """
        nodes = []
        for ID, node in self.graph.nodes.items():
            if ID == self.root.ID:
                nodes.append(root)
            else:
                nodes.append(node)
        Vsets = tolerantPartition(nodes, lambda u: u.x)
        Hsets = tolerantPartition(nodes, lambda u: u.y)
        return {'h': Hsets, 'v': Vsets}

    def getCCs(self, rs=None, ix=None):
        ccs = []
        for pc in self.pcs:
            ccs.extend(pc.buildCCs(rs=rs, ix=ix))
        return ccs

    def buildBufferNodesAndPCs(self, iel, dp, doBuildBufferNodes=True):
        """
        :param iel: ideal edge length for the graph
        :param dp: placement direction for this tree

        This method will be for the more complicated set of buffer
        nodes, which fit closer to the shape of the tree, if we decide
        to try that.
        See method by same name in TreePlacement class.

        We build buffer nodes to go on top of leaves and on the outside sides
        of nodes at the ends of ranks.

        The PCs will not only constrain the buffer nodes beside tree nodes, but
        will also maintain the shape of the tree: an alignment for each rank,
        sepcos maintaining both ordering and gaps within each rank, and either
        a rigid distribution on the ranks, or at least min gaps between them.
        """
        bns = []
        pcs = []
        if self.growthDir in Compass.vertical:
            axialDim, transDim = adg.YDIM, adg.XDIM
            axialCoord, transCoord = (lambda u: u.y), (lambda u: u.x)
            axialMeasure, transMeasure = (lambda u: u.h), (lambda u: u.w)
        else:
            axialDim, transDim = adg.XDIM, adg.YDIM
            axialCoord, transCoord = (lambda u: u.x), (lambda u: u.y)
            axialMeasure, transMeasure = (lambda u: u.w), (lambda u: u.h)
        # Build buffer nodes if requested.
        if doBuildBufferNodes:
            pad = iel/4.0
            def makeNode(x, X, y, Y, node=None):
                if node is not None:
                    u = node
                else:
                    u = Node()
                    u.ID = self.graph.getNextID()
                    u.fill = '#C0804080'
                    u.setIDAsLabel()
                u.w, u.h = X - x, Y - y
                u.x, u.y = x + u.w/2.0, y + u.h/2.0
                return u
            existingBNs = Deque(self.bufferNodes)
            # Pads on tops of leaves:
            for ID in self.leaves:
                leaf = self.nodes[ID]
                x, X, y, Y = leaf.boundingBoxxXyY()

                bn = existingBNs.popleft() if len(existingBNs) > 0 else None

                if self.growthDir == Compass.NORTH:
                    bn = makeNode(x, X, y - pad, y, node=bn)
                elif self.growthDir == Compass.SOUTH:
                    bn = makeNode(x, X, Y, Y + pad, node=bn)
                elif self.growthDir == Compass.EAST:
                    bn = makeNode(X, X + pad, y, Y, node=bn)
                else:
                    assert(self.growthDir == Compass.WEST)
                    bn = makeNode(x - pad, x, y, Y, node=bn)

                L, R = (leaf, bn) if self.growthDir in Compass.increasing else (bn, leaf)
                bns.append(bn)
                axialGap = pad/2.0 + axialMeasure(leaf)/2.0
                pcs.append(SepCo(axialDim, L, R, axialGap, exact=True))
                pcs.append(SepCo(transDim, L, R, 0, exact=True))
            # Pads on outsides of ranks:
            for i in range(1, self.depth):
                rank = self.nodesByRank[i]
                S = sorted(rank, key=transCoord)
                first, last = S[0], S[-1]
                x, X, y, Y = first.boundingBoxxXyY()
                u, U, v, V = last.boundingBoxxXyY()

                a = existingBNs.popleft() if len(existingBNs) > 0 else None
                b = existingBNs.popleft() if len(existingBNs) > 0 else None

                if self.growthDir in Compass.vertical:
                    a, b = makeNode(x - pad, x, y, Y, node=a), makeNode(U, U + pad, v, V, node=b)
                else:
                    assert(self.growthDir in Compass.horizontal)
                    a, b = makeNode(x, X, y - pad, y, node=a), makeNode(u, U, V, V + pad, node=b)

                bns.extend([a, b])
                firstGap = pad/2.0 + transMeasure(first)/2.0
                lastGap = pad/2.0 + transMeasure(last)/2.0
                pcs.append(SepCo(transDim, a, first, firstGap, exact=True))
                pcs.append(SepCo(axialDim, a, first, 0, exact=True))
                pcs.append(SepCo(transDim, last, b, lastGap, exact=True))
                pcs.append(SepCo(axialDim, last, b, 0, exact=True))
        # Generate the basic tree constraints.
        tallestNodes = []
        for i in range(self.depth):
            rank = self.nodesByRank[i]
            tallestNodes.append(max(rank, key=axialMeasure))
            # Align all nodes within the rank.
            pcs.append(AlignCo(axialDim, shapes=[(u, 0) for u in rank]))
            # Separate them.
            S = sorted(rank, key=transCoord)
            for L, R in zip(S[:-1], S[1:]):
                gap = transMeasure(L)/2.0 + iel/2.0 + transMeasure(R)/2.0
                pcs.append(SepCo(transDim, L, R, gap))
        # Rank separations:
        for i in range(self.depth - 1):
            A, B = tallestNodes[i:i+2]
            L, R = (A, B) if self.growthDir in Compass.increasing else (B, A)
            gap = axialMeasure(L)/2.0 + iel/2.0 + axialMeasure(R)/2.0
            pcs.append(SepCo(axialDim, L, R, gap, exact=self.holaConfig.RIGID_RANK_SEP))
        # Alignments with centre children:
        for p in self.nodes.values():
            if p == self.root and dp in Compass.cwOrds: continue
            Ch = p.getChildren()
            Ch.sort(key=transCoord)
            n = len(Ch)
            m = (n - (n % 2)) / 2
            if n % 2 == 1:
                pcs.append(AlignCo(transDim, shapes=[(p, 0), (Ch[m], 0)]))
            # Try flexible distributions on the "mirror triples?"
            if self.holaConfig.TRY_MIRROR_TRIPLES:
                for i in range(m):
                    a, z = Ch[i], Ch[n - 1 - i]
                    pcs.append(FlexDistCo(transDim, a, p, z))
        DEBUG = False
        if DEBUG:
            print
            print 'Tree constraints for tree rooted at %d:' % self.root.ID
            for pc in pcs:
                print pc
        self.pcs = pcs
        self.bufferNodes = bns
        return (bns, pcs)

    def writeFixedRelCo(self, root, extraNodes={}):
        """
        :param root: a Node
        :return: an adg.FixedRelativeConstraint, with the supplied root
                 substituted for this tree's root node
        """
        nodes = self.graph.nodes.copy()
        del nodes[self.root.ID]
        nodes[root.ID] = root
        for ID, u in extraNodes.items():
            nodes[ID] = u
        frc = FixedRelCo(nodes)
        return frc

    def insertIntoGraph(self, G):
        """
        :param G: a Graph
        :return: nothing

        Insert the tree graph into the given graph G, except for the root node.
        We assume the corresponding root node in G has the same ID.
        """
        H = self.graph
        for ID in H.nodes:
            if ID == self.root.ID: continue
            if G.nodes.has_key(ID):
                raise Exception("Inserting node %d. Graph already using that ID." % ID)
            G.nodes[ID] = H.nodes[ID]
        for edge in H.edges.values():
            s, t = edge.srcID, edge.tgtID
            G.addEdgeByNodeIDs(s, t)

    def size(self):
        return len(self.graph.nodes.keys())

    def buildOwnRoutingRig(self, rig, permissive=0):
        rig.addNodesAndEdges({self.root.ID:self.root}, {})
        self.addInternalNodesAndEdgesToRoutingRig(rig, permissive=permissive)
        return rig

    def addInternalNodesAndEdgesToRoutingRig(self, rig, core=None, permissive=0):
        # Get all nodes other than the root.
        nodes = {ID:node for ID, node in self.nodes.items() if ID != self.root.ID}
        # Get edges.
        edges = self.graph.edges
        # Set connection directions.
        connDirs = {}
        up = self.growthDir
        down = Compass.flip(up)
        upperDirs = Compass.libavoidVisibility[down]
        ordinalRootDirs = adg.ConnDirAll ^ upperDirs
        if permissive >= 2:
            lowerDirs = ordinalRootDirs
        else:
            lowerDirs = Compass.libavoidVisibility[up]
        dirs = (lowerDirs, upperDirs)
        rev = (upperDirs, lowerDirs)
        # If the tree is ordinally placed, then its root node
        # is in a different position than the root node belonging to the core.
        # So in order for the connectors in the router to get the proper endpts,
        # we temporarily move the root node to the position of the core root.
        coreroot = core.nodes[self.root.ID] if core is not None else self.root
        corerootpt = coreroot.centre()
        rootpt = self.root.centre()
        self.root.x, self.root.y = corerootpt
        # Compute just once the boolean to say whether we use permissive directions for
        # connection to the root node:
        permissiveRootDirections = (
            permissive >= 2 or
            (permissive == 1 and rootpt != corerootpt and len(self.nodesByRank[1]) == 1)
        )
        for rep, edge in edges.items():
            # Determine the ranks of the source and target ends, so that we can set the
            # conndirs in the right order.
            rs = self.rankByNodeID[edge.srcID]
            rt = self.rankByNodeID[edge.tgtID]
            # Use more permissive directions for root node if it is an ordinal placement.
            if permissiveRootDirections and self.root.ID in [edge.srcID, edge.tgtID]:
                d = (ordinalRootDirs, upperDirs) if rs < rt else (upperDirs, ordinalRootDirs)
            else:
                d = dirs if rs < rt else rev
            connDirs[rep] = d
        rig.addNodesAndEdges(nodes, edges, connDirs=connDirs)
        # Retore root position.
        self.root.x, self.root.y = rootpt

    def setEdgeRoutesInGraph(self, G):
        for e in self.graph.edges.values():
            f = G.getEdgeBtwNodes(e.src, e.tgt)
            #f = G.edges.get(rep, None)
            if f is not None:
                if f.srcID == e.srcID:
                    f.route = e.route
                else:
                    f.route = list(reversed(e.route))

    def getBoundingBoxXYWHWithoutRoot(self, growthDir=None):
        """
        Before calling this function the Tree should have been laid out,
        so the root node should be centred at (0, 0) and self.growthDir
        should be set.

        :param growthDir: a cardinal Compass direction
        :return: (ULCx, ULCy, w, h, u, v) where (ULCx, ULCy, w, h) gives the bounding
                 box after taking away the root node (but not taking away
                 inter-rank space between the root node and the next rank),
                 and where (u, v) is the vector from the centre of the root
                 node to the centre of the bounding box. Since the root is
                 at (0, 0), this is just the coords of the centre of the box.

        If you pass a growthDir argument then we use that. Otherwise we
        use self.growthDir.
        """
        assert(self.growthDir is not None)
        assert(self.root.x == 0 and self.root.y == 0)
        x, y, w, h = self.graph.boundingBoxXYWH()
        if self.growthDir == Compass.EAST:
            x += self.root.w
            w -= self.root.w
        elif self.growthDir == Compass.WEST:
            w -= self.root.w
        elif self.growthDir == Compass.SOUTH:
            y += self.root.h
            h -= self.root.h
        elif self.growthDir == Compass.NORTH:
            h -= self.root.h
        u, v = x + w/2.0, y + h/2.0
        if growthDir is not None:
            u, v = Compass.getRotationFunction(self.growthDir, growthDir)([u,v])
            if not Compass.sameDimension(self.growthDir, growthDir):
                # Dimensions of tree box swap.
                w, h = h, w
                # If root node is oblong, then tree box needs to be shifted
                # in the growth direction.
                if growthDir == Compass.EAST:
                    u += (self.root.w - self.root.h) / 2.0
                elif growthDir == Compass.WEST:
                    u -= (self.root.w - self.root.h) / 2.0
                elif growthDir == Compass.SOUTH:
                    v += (self.root.h - self.root.w) / 2.0
                else:
                    assert growthDir == Compass.NORTH
                    v -= (self.root.h - self.root.w) / 2.0
            # ULC is easily recomputed based on new centre.
            x, y = u - w/2.0, v - h/2.0
        return (x, y, w, h, u, v)

    def flip(self, growthDir):
        """
        :param growthDir: A compass direction, e.g. Compass.NORTH
        :return: nothing

        We filp the tree around the y-axis for NORTH and SOUTH growth,
        around the x-axis for EAST and WEST growth.
        """
        yAxis = growthDir in [Compass.NORTH, Compass.SOUTH]
        for node in self.nodes.values():
            if yAxis:
                node.x = -node.x
            else:
                node.y = -node.y
        self.lb, self.ub = (self.ub, self.lb)
        for k in self.boundsByRank:
            b = self.boundsByRank[k]
            self.boundsByRank[k] = [-b[1], -b[0]]

    def translate(self, vect, growthDir):
        """
        :param vect: a list or tuple of two numbers
        :return: nothing

        We add the vector to the coordinates of every node in the tree.
        """
        for node in self.nodes.values():
            node.x += vect[0]
            node.y += vect[1]
        c = 0 if growthDir in [Compass.NORTH, Compass.SOUTH] else 1
        self.lb += vect[c]
        self.ub += vect[c]
        for k in self.boundsByRank:
            b = self.boundsByRank[k]
            self.boundsByRank[k] = [b[0] + vect[c], b[1] + vect[c]]

    def rotateTo(self, growthDir):
        """
        :param growthDir: a cardinal Compass direction
        :return: nothing

        We rotate the tree around (0,0) so that it grows in the given direction.
        We assume that by now the tree's self.growthDir has been defined.
        If not, we raise an exception!
        """
        r = Compass.getRotationFunction(self.growthDir, growthDir)
        for node in self.nodes.values():
            node.x, node.y = r((node.x, node.y))
        self.growthDir = growthDir

    def getBounds(self, rank, nodeSep):
        if self.boundaryInfinite:
            if self.boundaryTight:
                bds = self.boundsByRank.get(rank, [0,0])
            else:
                bds = [self.lb, self.ub]
        else:
            bds = self.boundsByRank.get(rank, None)
            if not self.boundaryTight and bds is not None:
                bds = [self.lb, self.ub]
        h = nodeSep
        if bds is not None:
            bds = [bds[0] - h, bds[1] + h]
        return bds

    def setRank(self, node, rank):
        self.rankByNodeID[node.ID] = rank
        R = self.nodesByRank.get(rank, [])
        R.append(node)
        self.nodesByRank[rank] = R

    def getCTrees(self):
        """
        The "C-Trees" are the connected components of the graph obtained
        by deleting the root node of this tree.

        But for functional C-Trees we don't actually need to delete any
        nodes, we just need to make sure the new Trees have the right root
        node, depth, breadth, and mapping from rank to nodes in that rank.
        The set of leaves is the same, since being in the subtree and
        being a leaf of the original tree are together equivalent to
        being a leaf of the subtree.
        """
        cTrees = []
        children = self.root.getChildren()
        for c in children:
            t = Tree(self.graph, c)
            t.setHolaConfig(self.holaConfig)
            cTrees.append(t)
        return cTrees

    def leavesOfRank(self, r):
        L = []
        N = self.nodesByRank[r]
        for node in N:
            if node.ID in self.leaves:
                L.append(node)
        return L

    def nonleavesOfRank(self, r):
        L = []
        N = self.nodesByRank[r]
        for node in N:
            if node.ID not in self.leaves:
                L.append(node)
        return L

    def computeIsomString(self):
        """
        Compute a string which uniquely represents the
        isomorphism class of this tree.
        This is the core idea from Manning & Atallah 1985.
        """
        isomNumber = {}
        isomTuple = {}
        isomTupleString = {}
        levelIsomStrings = []
        # Assign isomNumber 0 to all leaves.
        for ID in self.leaves: isomNumber[ID] = 0
        # Let L be a list of the deepest leaves.
        d = self.depth
        L = self.leavesOfRank(d-1)
        # Compute the isomstring for each level.
        for r in range(d-2, -1, -1):
            N = self.nonleavesOfRank(r)
            for node in N: isomTuple[node.ID] = []
            # For each leaf v in L, insert its number into
            # the tuple of its parent node.
            for v in L:
                p = self.parent[v.ID]
                n = isomNumber[v.ID]
                isomTuple[p.ID].append(n)
            # Create strings for isomTuples.
            for u in N:
                isomTupleString[u.ID] = ','.join([
                    str(n) for n in sorted(isomTuple[u.ID])
                ])
            # Sort the nonleaves of rank r by tuple.
            N.sort(key=lambda u: isomTupleString[u.ID])
            # Now can write isom string for this level.
            levelIsomStrings.append( ';'.join(
                isomTupleString[u.ID] for u in N
            ))
            # Compute next L.
            L = []
            nodesByTuple = {}
            for u in N:
                t = isomTupleString[u.ID]
                A = nodesByTuple.get(t,[])
                A.append(u)
                nodesByTuple[t] = A
            distinctTuples = nodesByTuple.keys()
            distinctTuples.sort()
            for k in range(len(distinctTuples)):
                t = distinctTuples[k]
                nodes = nodesByTuple[t]
                for u in nodes:
                    isomNumber[u.ID] = k+1
                    L.append(u)
            leaves = self.leavesOfRank(r)
            leaves.extend(L)
            L = leaves
        # Join the level strings for the full tree string.
        return ':'.join(levelIsomStrings)

    def computeOrderingNumbers(self):
        """
        Assign each node an "ordering number," storing these in self.orderingNumbers.
        These numbers are to be used by the self.cmpRankmates function in order to sort
        the nodes of each rank when it comes time to generate constraints for the tree.

        :return: nothing
        """
        # Leaves are simple.
        if self.depth == 1:
            self.isSymmetric = True
            return
        # Proceed for nonleaves:
        C = self.getCTrees()
        # Layout the C-Trees, recursively.
        for t in C: t.computeOrderingNumbers()
        # Sort the C-trees into isomorphism classes.
        classes = {}
        for t in C:
            isomstr = t.computeIsomString()
            A = classes.get(isomstr, [])
            A.append(t)
            classes[isomstr] = A
        # Now sort the classes.
        isoms = classes.keys()
        def isomCmp(I,J):
            c = classes[I]; d = classes[J]
            cd = c[0].depth; dd = d[0].depth
            cb = c[0].breadth; db = d[0].breadth
            # Put narrower trees first.
            if cb > db: return 1
            if cb < db: return -1
            # For same breadth, put shallower trees first.
            if cd > dd: return 1
            if cd < dd: return -1
            # Otherwise just compare the isomorphism strings for some
            # way to make this relation deterministic.
            if I < J: return -1
            if I > J: return 1
            return 0
        isoms.sort(cmp=isomCmp)
        # Which classes have odd order?
        oddOrder = {}
        for I in isoms:
            c = classes[I]
            if len(c) % 2 == 1:
                oddOrder[I] = c
        # Determine whether our layout is going to be symmetric or not.
        numOdd = len(oddOrder)
        haveCentralTree = False
        # If there are no odd-order classes, then we are symmetric.
        if numOdd == 0:
            self.isSymmetric = True
        # If there are two or more odd-order classes, then we are not symmetric.
        elif numOdd > 1:
            self.isSymmetric = False
        # Else there is exactly one odd-order class.
        # In this case we are symmetric if and only if (any representative of) the one
        # odd order class is symmetric.
        else:
            self.isSymmetric = oddOrder.values()[0][0].isSymmetric
            # For symmetric layout the trees of odd-order class need to go in the centre,
            # so we put them first in the list, since we work our way outward from the centre
            # when placing the trees.
            oddIsom = oddOrder.keys()[0]
            isoms.remove(oddIsom)
            isoms.insert(0, oddIsom)
        # Now order the trees alternating around the centre, flipping the trees that get placed
        # on the left hand side.
        signedtrees = Deque()
        nextsign = 1
        for I in isoms:
            c = classes[I]
            for t in c:
                op = signedtrees.append if nextsign == 1 else signedtrees.appendleft
                op((nextsign, t))
                nextsign *= -1
        # Set ordering numbers.
        for i, st in enumerate(signedtrees):
            s, t = st
            self.orderingNumbers[t.root.ID] = i
            for ID, j in t.orderingNumbers.items():
                self.orderingNumbers[ID] = s*j

    def writeConstraints(self, iel, dp, rigidRankSep=False):
        pcs = []
        if self.growthDir in Compass.vertical:
            axialDim, transDim = adg.YDIM, adg.XDIM
            axialCoord, transCoord = (lambda u: u.y), (lambda u: u.x)
            axialMeasure, transMeasure = (lambda u: u.h), (lambda u: u.w)
        else:
            axialDim, transDim = adg.XDIM, adg.YDIM
            axialCoord, transCoord = (lambda u: u.x), (lambda u: u.y)
            axialMeasure, transMeasure = (lambda u: u.w), (lambda u: u.h)
        doReverse = self.growthDir in Compass.decreasing
        # Generate the basic tree constraints.
        tallestNodes = []
        for i in range(1, self.depth):
            rank = self.nodesByRank[i]
            tallestNodes.append(max(rank, key=axialMeasure))
            # Align all nodes within the rank.
            pcs.append(AlignCo(axialDim, shapes=[(u, 0) for u in rank]))
            # Separate them.
            S = sorted(rank, cmp=self.cmpRankmates, reverse=doReverse)
            for L, R in zip(S[:-1], S[1:]):
                gap = transMeasure(L)/2.0 + iel/2.0 + transMeasure(R)/2.0
                pcs.append(SepCo(transDim, L, R, gap))
        # Rank separations:
        for i in range(self.depth - 1):
            A, B = tallestNodes[i:i+2]
            L, R = (A, B) if self.growthDir in Compass.increasing else (B, A)
            gap = axialMeasure(L)/2.0 + iel/2.0 + axialMeasure(R)/2.0
            pcs.append(SepCo(axialDim, L, R, gap, exact=rigidRankSep))
        # Alignments with centre children:
        for p in self.nodes.values():
            if p == self.root and dp in Compass.cwOrds: continue
            Ch = p.getChildren()
            Ch.sort(cmp=self.cmpRankmates, reverse=doReverse)
            n = len(Ch)
            m = (n - (n % 2)) / 2
            if n % 2 == 1:
                pcs.append(AlignCo(transDim, shapes=[(p, 0), (Ch[m], 0)]))
            # Try flexible distributions on the "mirror triples?"
            if self.holaConfig.TRY_MIRROR_TRIPLES:
                for i in range(m):
                    a, z = Ch[i], Ch[n - 1 - i]
                    pcs.append(DistCo(transDim, [a, p, z],
                                      iel*2*(m-i),  # actually shouldn't matter
                                      flexible=True))
        DEBUG = False
        if DEBUG:
            print
            print 'Tree constraints for tree rooted at %d:' % self.root.ID
            for pc in pcs:
                print pc
        self.pcs = pcs
        return pcs

    def symmetricLayout(self, growthDir, nodeSep, rankSep):
        """
        Alternative to the 'buildSymmetricTreeConstraints' method.

        Instead of computing symmetric tree constraints, you can just compute a symmetric
        layout.

        :growthDir: a compass direction, being a field of the Compass enum,
                    e.g. Compass.NORTH
        :nodeSep: minimal gap between nodes on the same rank
        :rankSep: separation between ranks
        :return: nothing
        """
        # Save the growth direction.
        self.growthDir = growthDir
        # Initialise root position to zero.
        self.root.x = 0
        self.root.y = 0
        # Initialise rank bounds to zero.
        for k in range(self.depth):
            self.boundsByRank[k] = [0, 0]
        # Initialise rank 0 bounds for root node.
        half = self.root.w/2.0 if growthDir in [Compass.NORTH, Compass.SOUTH] else self.root.h/2.0
        self.lb = -half
        self.ub = half
        self.boundsByRank[0] = [self.lb, self.ub]
        # Leaves are simple.
        if self.depth == 1:
            self.isSymmetric = True
            return
        # Proceed for nonleaves:
        C = self.getCTrees()
        # Layout the C-Trees, recursively.
        for t in C: t.symmetricLayout(growthDir, nodeSep, rankSep)
        # Sort the C-trees into isomorphism classes.
        classes = {}
        for t in C:
            isomstr = t.computeIsomString()
            A = classes.get(isomstr, [])
            A.append(t)
            classes[isomstr] = A
        # Now sort the classes.
        isoms = classes.keys()
        def isomCmp(I,J):
            c = classes[I]; d = classes[J]
            cd = c[0].depth; dd = d[0].depth
            cb = c[0].breadth; db = d[0].breadth
            if cb > db: return 1
            if cb < db: return -1
            if cd > dd: return 1
            if cd < dd: return -1
            # Otherwise just compare the isomorphism strings for some
            # way to make this relation deterministic.
            if I < J: return -1
            if I > J: return 1
            return 0
        isoms.sort(cmp=isomCmp)
        # Which classes have odd order?
        oddOrder = {}
        for I in isoms:
            c = classes[I]
            if len(c) % 2 == 1:
                oddOrder[I] = c
        # Determine whether our layout is going to be symmetric or not.
        numOdd = len(oddOrder)
        haveCentralTree = False
        # If there are no odd-order classes, then we are symmetric.
        if numOdd == 0:
            self.isSymmetric = True
        # If there are two or more odd-order classes, then we are not symmetric.
        elif numOdd > 1:
            self.isSymmetric = False
        # Else there is exactly one odd-order class.
        # In this case we are symmetric if and only if (any representative of) the one
        # odd order class is symmetric.
        else:
            self.isSymmetric = oddOrder.values()[0][0].isSymmetric
            # For symmetric layout the trees of odd-order class need to go in the centre,
            # so we put them first in the list, since we work our way outward from the centre
            # when placing the trees.
            haveCentralTree = True
            oddIsom = oddOrder.keys()[0]
            isoms.remove(oddIsom)
            isoms.insert(0, oddIsom)
        # Now place the c-trees alternating around the centre.
        # For this operation this tree must have "infinite boundary" at the
        # centre line. But we preserve the configured value and restore it
        # when we are done.
        givenBoundaryInfiniteValue = self.boundaryInfinite
        self.boundaryInfinite = True
        positiveNext = True
        mustPlaceCentralTree = haveCentralTree
        baseTrans = {
            Compass.NORTH: [0, -rankSep],
            Compass.EAST: [rankSep, 0],
            Compass.SOUTH: [0, rankSep],
            Compass.WEST: [-rankSep, 0]
        }[growthDir]
        ns = nodeSep
        for I in isoms:
            c = classes[I]
            for t in c:
                if mustPlaceCentralTree:
                    t.translate(baseTrans, growthDir)
                    lb, ub = [0, 0]
                    for i in range(self.depth):
                        if 1 <= i and i <= t.depth:
                            tlb, tub = t.boundsByRank[i-1]
                            self.boundsByRank[i] = [tlb, tub]
                            lb = min(lb, tlb)
                            ub = max(ub, tub)
                    self.lb
                    self.ub
                    mustPlaceCentralTree = False
                else:
                    # Set up based on whether we're on the positive or
                    # the negative side.
                    a = 1 if positiveNext else 0
                    b = 0 if positiveNext else 1
                    rootPosChooser = max if positiveNext else min
                    if not positiveNext: t.flip(growthDir)
                    positiveNext = not positiveNext
                    # Compute the position where each rank would like the
                    # root to go.
                    inf = t.boundaryInfinite
                    t.boundaryInfinite = True
                    rootPosPerRank = [
                        self.getBounds(i+1, ns)[a] - t.getBounds(i, ns)[b]
                        for i in range(t.depth)
                    ]
                    t.boundaryInfinite = inf
                    # Take max or min, according to which side we are on.
                    rootPos = rootPosChooser(rootPosPerRank)
                    # Compute the translation for the tree.
                    if growthDir in [Compass.NORTH, Compass.SOUTH]:
                        trans = [rootPos, baseTrans[1]]
                    else:
                        trans = [baseTrans[0], rootPos]
                    # Translate it.
                    t.translate(trans, growthDir)
                    # Update upper bounds if a == 1, lower if a == 0.
                    extreme = 0
                    for i in range(self.depth):
                        if 1 <= i and i <= t.depth:
                            bd = t.boundsByRank[i-1][a]
                            self.boundsByRank[i][a] = bd
                        else:
                            bd = self.boundsByRank[i][a]
                        extreme = rootPosChooser(extreme, bd)
                    if a == 1:
                        self.ub = bd
                    else:
                        self.lb = bd
        # Restore configuration
        self.boundaryInfinite = givenBoundaryInfiniteValue
