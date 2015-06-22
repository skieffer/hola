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

from graphs import Graph
from graphs import Edge

class Stem:

    def __init__(self,leaf,root):
        self.leaf = leaf
        self.root = root

    def addSelfToGraph(self, H):
        """
        First make sure the graph H has nodes of IDs matching those of
        the endpts of this stem, creating them from "excised copies" if
        necessary. Set the treeSerialNo's of the nodes, and connect
        them in H.
        """
        leaf = H.getNode(self.leaf.ID)
        if not leaf:
            leaf = self.leaf.excisedCopy()
            leaf.treeSerialNo = 0
            H.addNode(leaf)
        root = H.getNode(self.root.ID)
        if not root:
            root = self.root.excisedCopy()
            H.addNode(root)
        # Must always update root node's serial number to be the largest.
        root.treeSerialNo = H.takeNextTreeSerialNo()
        # We always make the root end the source, so that the directed edges
        # in the tree will flow from root to leaves.
        E = Edge(srcID=root.ID, tgtID=leaf.ID)
        H.addEdge(E)


def stemsFromLeaves(L):
    """
    L: {ID:Node} dict
    return: list of Stems
    """
    stems = []
    for ID in L:
        leaf = L[ID]
        el = leaf.edgeList()
        root = el[0].otherEnd(leaf)
        s = Stem(leaf, root)
        stems.append(s)
    return stems


class NodeBuckets:
    """
    For sorting all the nodes of a graph into "buckets" according
    to their degree.
    You can take all the leaves (degree-1 nodes) and you can move
    a node from one bucket to another.
    """

    def __init__(self,G):
        self.graph = G
        self.buckets = []
        self.M = G.maxDeg
        for d in range(0,self.M+1):
            self.buckets.append({})
        nodes = G.nodes
        for ID in nodes:
            node = nodes[ID]
            deg = node.degree
            self.buckets[deg][ID] = node

    def takeLeaves(self):
        L = self.buckets[1]
        self.buckets[1] = {}
        return L

    def hasLeaves(self):
        return len(self.buckets[1].keys()) > 0

    def cutOneStem(self):
        """
        If there are any leaves, make a Stem object to represent one of
        them, cut the leaf from the graph, and adjust buckets.
        """
        stem = None
        if self.hasLeaves():
            # Choose a leaf.
            k0 = self.buckets[1].keys()[0]
            leaf = self.buckets[1][k0]
            # Make a stem for it.
            el = leaf.edgeList()
            edge = el[0]
            root = edge.otherEnd(leaf)
            stem = Stem(leaf, root)
            # Move nodes up in the bucket stack.
            self.moveNode(leaf.ID, 1, 0)
            rootDeg = root.degree
            self.moveNode(root.ID, rootDeg, rootDeg - 1)
            # Cut from graph.
            self.graph.severEdge(edge)
            self.graph.deleteNode(leaf)
        return stem

    def moveNode(self,ID,oldDeg,newDeg):
        """
        Move node of given ID from old degree to new degree.
        Fails quietly if there is no node of given ID in
        oldDeg bucket.
        """
        try:
            node = self.buckets[oldDeg][ID]
            self.buckets[newDeg][ID] = node
            del self.buckets[oldDeg][ID]
        except: pass


def prune(G):
    """
    G: an instance of the Graph class
    return: a list of Graphs in which the first element is the trunk
            and the remaining elements are the trees.

    The exception is that if G was already a tree then we return
    it as the first and only component of the list.

    The treeSerialNo fields of the nodes in all non-trunk components will have
    been set, and the identifyRootNode method called on each, so their
    self.rootNode fields are already set to the root node of the tree.

    Note: it IS necessary to take all the nodes which are leaves at a given moment
    all together, and then take the next batch of leaves.
    If you just pick any leaf, one by one, then the root can wind up
    anywhere in the entire tree, and you don't actually find the "centre" node.
    """
    B = NodeBuckets(G)
    H = Graph()
    L = B.takeLeaves()
    while L:
        S = stemsFromLeaves(L)
        G.severNodes(L,buckets=B)
        if G.isEmpty():
            # This case only arises when G was a tree with a double centre
            # U--V. At this point we must have:
            #     S == [Stem(U,V), Stem(V,U)]
            # and we must only add one (doesn't matter which) of those stems
            # to H.
            S = S[:1]
        for stem in S: stem.addSelfToGraph(H)
        L = B.takeLeaves()
    C = H.connComps()
    for c in C:
        c.identifyRootNode()
    # G will now have 3 or more nodes iff it was
    # not in fact a tree, to begin with.
    if G.getNumNodes() >= 3:
        C.insert(0,G)
        # It will be helpful to let each root node in G know that it is a root node.
        for c in C[1:]:
            rp = c.rootNode
            r = G.nodes[rp.ID]
            r.isRootNode = True
    return C

