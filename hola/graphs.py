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

import sys
import random
import json
import re
import copy
import math
from collections import deque as Deque

import adaptagrams.adaptagrams as adg
from logging import LogLevel
from util import NearbyObjectFinder
from constraints import *
import svg


class idDispenser:

    def __init__(self, firstID):
        self.nextID = firstID

    def takeNext(self):
        n = self.nextID
        self.nextID = n + 1
        return n

    def noteIDInUse(self, ID):
        self.nextID = max(self.nextID, ID + 1)

    def reset(self):
        self.nextID = 0

ID_DISP = idDispenser(0)

def buildGraphFromTGF(tgf, w=lambda a: 30, h=lambda a: 30):
    """
    :param tgf: graph defined in Trivial Graph Format
    :param w: function that returns a width for each node ID
    :param h: function that returns a height for each node ID
    :return: a Graph object
    """
    G = Graph()
    lines = tgf.split('\n')
    mode = 0
    def intIfInt(s):
        try:
            s = int(s)
        except:
            pass
        return s
    for line in lines:
        if len(line) == 0: continue
        if line.strip() == "#":
            mode = 1
        elif mode == 0:
            ID = line.strip()
            ID = intIfInt(ID)
            u = Node()
            u.ID = ID
            u.w = w(ID)
            u.h = h(ID)
            G.addNode(u)
        elif mode == 1:
            srcID, tgtID = [intIfInt(s) for s in line.strip().split()]
            e = Edge(srcID, tgtID)
            G.addEdge(e)
    return G

class Graph:

    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.nodeConf = NodeConfig(self)
        self.extraPCs = []
        self.PCHolders = []
        self.maxDeg = 0
        self.nextTreeSerialNo = 1
        self.iel = 0

        # For debugging:
        self.logger = None
        self.applyProjSeqDebugNumber = 0

        self.nodePosStack = []

        # Dictionary for storing shortest path lengths (in graph-theoretic hops), to be looked up
        # by node IDs. I.e. distance between nodes u and v will be self.shortestPaths[u.ID][v.ID].
        self.shortestPaths = {}
        # Likewise, but for the "G" connectivity matrix:
        self.connectionMatrix = {}
        # Boolean to control how we compute stress:
        self.COMPUTE_STRESS_LOCAL = False

    EDGENODE_THICKNESS = 10

    def addPaddingToNodes(self, xPad, yPad):
        # padding may be positive or negative
        for node in self.nodes.values():
            node.addPadding(xPad, yPad)

    def randomLayout(self):
        d = self.computeAvgNodeDim()
        n = len(self.nodes)
        R = 3*d*math.sqrt(n)
        for node in self.nodes.values():
            node.x = random.uniform(0, R)
            node.y = random.uniform(0, R)

    def computeShortestPaths(self, iel):
        """
        Compute the shortests paths between all pairs of nodes, and
        store the result in self.shortestPaths. Must be recomputed if
        the graph structure is changed (i.e. nodes or edges are added
        or taken away).

        We also compute the connection matrix while we're at it.
        :return: nothing
        """
        rs, es, ix = self.writeRsEsIx()
        # Prepare inverse mapping to ix.
        ixinv = {}
        for node in self.nodes.values():
            ixinv[ix(node.ID)] = node.ID
        #
        op = False
        alg = adg.ConstrainedFDLayout(rs, es, iel, op)
        d = alg.readLinearD()
        g = alg.readLinearG()
        n = len(rs)
        sp = {}
        cm = {}
        for i in range(n):
            spRow = {}
            cmRow = {}
            for j in range(n):
                spRow[ixinv[j]] = d[n*i + j]
                cmRow[ixinv[j]] = g[n*i + j]
            sp[ixinv[i]] = spRow
            cm[ixinv[i]] = cmRow
        self.shortestPaths = sp
        self.connectionMatrix = cm
        self.COMPUTE_STRESS_LOCAL = True

    def computeStress(self, iel, logger=None):
        if self.COMPUTE_STRESS_LOCAL:
            return self.computeStressLocal(logger=logger)
        else:
            return self.computeStressOutsourced(iel, logger=logger)

    def computeStressLocal(self, logger=None):
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('compute_stress_local')
        stress = 0
        nodes = self.nodes.values()
        for i in range(len(nodes)-1):
            u = nodes[i]
            for j in range(i+1, len(nodes)):
                v = nodes[j]
                # If any nodes have been added since the shortest paths matrix was computed,
                # we assume those nodes are not important, or else you should not be computing
                # the stress now. So we just skip those cases.
                p = self.connectionMatrix.get(u.ID, {}).get(v.ID, None)
                if p is None: continue
                #p = self.connectionMatrix[u.ID][v.ID]
                if p == 0: continue  # no forces btw disconnected parts of the graph
                rx, ry = u.x - v.x, u.y - v.y
                l = math.sqrt(rx*rx + ry*ry)
                d = self.shortestPaths[u.ID][v.ID]
                if l > d and p > 1: continue  # no attractive forces required
                d2 = d*d
                rl = d - l
                s = rl*rl/float(d2)
                stress += s
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        return stress

    def computeStressOutsourced(self, iel, logger=None):
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('setup_FD_layout')
        alg = self.setupFDLayout(iel, logger=logger)[0]
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
            logger.startNewTimer('adg_compute_stress')
        stress = alg.computeStress()
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        return stress

    def writeShortestPathsMatrix(self):
        s = ''
        IDs = sorted([u.ID for u in self.nodes.values()])
        n = len(IDs)
        M = max([len(str(ID)) for ID in IDs])
        for i in IDs:
            for j in IDs:
                M = max(M, len(str(self.shortestPaths[i][j])))
        fmt = '%%%dd'%M
        headingFmt = ' '*M + '  ' + ' '.join([fmt]*n) + '\n'
        heading = headingFmt % tuple(IDs)
        s += heading
        s += '-'*len(heading) + '\n'
        rowFmt =  fmt + "| " + ' '.join([fmt]*n) + '\n'
        for ID in IDs:
            s += rowFmt % (
                (ID,) + tuple([self.shortestPaths[ID][j] for j in IDs])
            )
        return s

    def setMinSize(self, w, h):
        for node in self.nodes.values():
            node.w = max(w, node.w)
            node.h = max(h, node.h)

    def addPCs(self, pcs):
        """
        :param pcs: list of python constraints
        :return: nothing

        We add the pcs to self.extraPCs.
        """
        self.extraPCs.extend(pcs)

    def addPCHolders(self, pchs):
        self.PCHolders.extend(pchs)

    def getAdjacencyMatrix(self):
        """
        :return: a sparse adjacency matrix am

        For IDs u, v of adjacent nodes, both am[u][v] and am[v][u] return True.
        But if the nodes are not adjacent, then am[u][v] is NOT DEFINED.
        However, am[u] IS DEFINED for all node IDs u.
        """
        am = {}
        for ID in self.nodes.keys():
            am[ID] = {}
        for edge in self.edges.values():
            am[edge.srcID][edge.tgtID] = True
            am[edge.tgtID][edge.srcID] = True
        return am

    def buildRigForEdgesNeedingRoute(self):
        rig = RoutingRig({
            'routing': adg.OrthogonalRouting,
            'nudgeOrthogonalSegmentsConnectedToShapes': True,
            'nudgeSharedPathsWithCommonEndPoint': True
        })
        edges = {r:e for r, e in self.edges.items() if e.route is None}
        rig.addNodesAndEdges(self.nodes, edges)
        return rig

    def doUseYEdColour(self):
        for node in self.nodes.values():
            node.useYEdColour = True

    def eraseAllLabels(self):
        for node in self.nodes.values():
            node.label = ""

    def clearAllRoutes(self):
        for edge in self.edges.values():
            edge.clearRoute()

    def buildRoutes(self):
        for edge in self.edges.values():
            edge.buildRoute()

    def setPosesInCorrespNodes(self, H):
        """
        :param H: another Graph object
        :return: nothing

        For each of own nodes, see if H has one of the same ID. If so, set the position
        of the node in H to equal the position of the node in self.
        """
        for ID in self.nodes.keys():
            h = H.nodes.get(ID, None)
            if h is not None:
                g = self.nodes[ID]
                h.x, h.y = g.x, g.y

    def setRoutesInCorrespEdges(self, H):
        """
        :param H: another Graph object
        :return: nothing

        For each of own edges, see if H has one of the same description.
        If so, set the route of the edge in H to equal the route of the edge in self.
        """
        for dsc in self.edges.keys():
            h = H.edges.get(dsc, None)
            if h is not None:
                g = self.edges[dsc]
                h.route = g.route
                h.routePoints = g.routePoints
                h.bends = g.bends

    def setLogger(self, logger):
        self.logger = logger
        logger.graph = self

    def __repr__(self):
        s = ''
        s += 'Graph:\n'
        s += '    Nodes:\n'
        for ID in sorted(self.nodes.keys()):
            N = self.nodes[ID]
            s += ' '*8 + '%2d: %s\n'%(ID, repr(N))
        s += '    Edges:\n'
        for rep in sorted(self.edges.keys()):
            s += ' '*8 + rep + '\n'
        return s

    def copy(self):
        H = Graph()
        for ID in self.nodes: H.nodes[ID] = self.nodes[ID]
        for rep in self.edges: H.edges[rep] = self.edges[rep]
        H.extraPCs = self.extraPCs[:]
        H.PCHolders = self.PCHolders[:]
        H.nodeConf = self.nodeConf.copy(H)
        H.maxDeg = self.maxDeg
        return H

    def translateBy(self, vect):
        dx, dy = vect
        for node in self.nodes.values():
            node.x += dx
            node.y += dy
        for edge in self.edges.values():
            edge.translateBy(vect)

    def rotate90cw(self, iel_for_constraints, iel_for_shake, logger=None, nbr=False):
        "Since nodes may be oblong rects, need to rotate constraints too, and shake."
        for node in self.nodes.values():
            node.x, node.y = -node.y, node.x
        for edge in self.edges.values():
            edge.rotate90cw()
        self.nodeConf.rotate90cw()
        for pc in self.extraPCs:
            pc.rotateCW(1)
        for pch in self.PCHolders:
            pch.rotatePCsCW(1, iel_for_constraints)
        if logger is not None and logger.level >= LogLevel.STAGE_GRAPHS:
            logger.writeGML("_13_0_rotateAndTranslate_before_shake", graph=self)
        self.shakeWithSolidEdges(iel_for_shake, useNeighbourStress=nbr)

    def rotate90acw(self, iel_for_constraints, iel_for_shake, logger=None, nbr=False):
        "Since nodes may be oblong rects, need to rotate constraints too, and shake."
        for node in self.nodes.values():
            node.x, node.y = node.y, -node.x
        for edge in self.edges.values():
            edge.rotate90acw()
        self.nodeConf.rotate90acw()
        for pc in self.extraPCs:
            pc.rotateCW(3)
        for pch in self.PCHolders:
            pch.rotatePCsCW(3, iel_for_constraints)
        if logger is not None and logger.level >= LogLevel.STAGE_GRAPHS:
            logger.writeGML("_13_0_rotateAndTranslate_before_shake", graph=self)
        self.shakeWithSolidEdges(iel_for_shake, useNeighbourStress=nbr)

    def rotate180(self, iel_for_constraints, iel_for_shake):
        for node in self.nodes.values():
            node.x, node.y = -node.x, -node.y
        for edge in self.edges.values():
            edge.rotate180()
        self.nodeConf.rotate180()
        for pc in self.extraPCs:
            pc.rotateCW(2)
        for pch in self.PCHolders:
            pch.rotatePCsCW(2, iel_for_constraints)

    def zeroAllPositions(self):
        for node in self.nodes.values():
            node.x, node.y = (0, 0)

    def boundingBoxxXyY(self, ignore=[], includeBends=False):
        """
        :return: bounding box in the form (x, X, y, Y) giving extreme coords
        """
        nodes = self.nodes.values()
        nodes = filter(lambda u: u not in ignore, nodes)
        bbox = nodes[0].boundingBoxxXyY()
        for node in nodes[1:]:
            u, U, v, V = node.boundingBoxxXyY()
            x, X, y, Y = bbox
            bbox = (min(x, u), max(X, U), min(y, v), max(Y, V))
        if includeBends:
            for edge in self.edges.values():
                u, U, v, V = edge.boundingBoxxXyY()
                x, X, y, Y = bbox
                bbox = (min(x, u), max(X, U), min(y, v), max(Y, V))
        return bbox

    def boundingBoxXYWH(self, ignore=[], includeBends=False):
        """
        :return: bounding box of the graph as (ULCx, ULCy, W, H)
        """
        x, X, y, Y = self.boundingBoxxXyY(ignore=ignore, includeBends=includeBends)
        return (x, y, X - x, Y - y)

    def bboxPerimeter(self):
        x, y, w, h = self.boundingBoxXYWH()
        return 2*(w + h)

    def buildBendPoints(self):
        """
        Look through the edges, and check their 'route' attribute.
        If the route is defined and length 3 or greater, create a new Node
        object for each internal point on the route.

        As a(n important) side effect, the sequence of bend point nodes
        created for each edge are stashed, in order, in the edge as
        self.routePoints.

        NB: One complication is that if there is a single point (x0, y0) at
        which two or more edges each have a bendpoint, then we are only
        going to create /one/ node at that point, and /all/ the edges are
        going to have it as one of their routePoints.

        :return: dict of the bend point nodes
        """
        bendpoints = {}

        bpsByXY = NearbyObjectFinder(0.5)
        for edge in self.edges.values():
            route = edge.route
            if route is not None and len(route) >= 3:
                routePoints = []
                for pt in route[1:-1]:
                    x, y = pt
                    bp = bpsByXY.findObject(x, y)
                    if bp is None:
                        ID = ID_DISP.takeNext()
                        bp = newBendNode(ID, x, y)
                        bp.setIDAsLabel()
                        bendpoints[bp.ID] = bp
                        bpsByXY.addObject(x, y, bp)
                    routePoints.append(bp)
                edge.routePoints = routePoints
        return bendpoints

    def setIDsAsLabels(self):
        for node in self.nodes.values():
            node.setIDAsLabel()

    def setIndicesAsLabels(self, ix):
        for node in self.nodes.values():
            node.label = str(ix(node))

    def getExtremeIDs(self):
        IDs = self.nodes.keys()
        minID = min(IDs)
        maxID = max(IDs)
        return (minID, maxID)

    def registerMaxID(self):
        """
        Tell the global ID dispenser about the max ID being used in this
        graph.
        :return: nothing
        """
        m, M = self.getExtremeIDs()
        ID_DISP.noteIDInUse(M)

    def getNextID(self):
        #minID, maxID = self.getExtremeIDs()
        #nextID = maxID + 1
        #return nextID
        return ID_DISP.takeNext()

    def getIDdispenser(self):
        #nextID = self.getNextID()
        #idd = idDispenser(nextID)
        #return idd
        return ID_DISP

    def createAndAddNewBendNode(self, x, y):
        # Add a bend node.
        ID = self.getNextID()
        bp = newBendNode(ID, x, y)
        self.addNode(bp)
        return bp

    def createAndAddNewEdge(self, srcID, tgtID):
        e = Edge(srcID, tgtID)
        self.addEdge(e)
        return e

    def shiftIDs(self, shift):
        """
        Add a constant to the IDs of all nodes.
        Adjust the edges accordingly.
        """
        nodes = {}
        edges = {}
        for ID in self.nodes:
            N = self.nodes[ID]
            N.ID += shift
            nodes[N.ID] = N
        for rep in self.edges:
            E = self.edges[rep]
            E.srcID += shift
            E.tgtID += shift
            edges[repr(E)] = E
        self.nodes = nodes
        self.edges = edges

    def takeNextTreeSerialNo(self):
        n = self.nextTreeSerialNo
        self.nextTreeSerialNo = n + 1
        return n

    def initFromGML(self, d):
        self.directed = (d.get('directed', None) == 1)
        NL = d['nodeList']
        for N in NL:
            self.addNode(N)
        EL = d['edgeList']
        for E in EL:
            self.addEdge(E) 

    def addNode(self, N):
        N.setGraph(self)
        self.nodes[N.ID] = N

    def addEdge(self, E):
        # We assume all nodes have been added before you
        # try to add any edges.
        E.setGraph(self)
        self.edges[repr(E)] = E
        self.maxDeg = max(self.maxDeg, E.src.degree, E.tgt.degree)

    def addEdgeByNodeIDs(self, srcID, tgtID):
        self.addEdge(Edge(srcID=srcID, tgtID=tgtID))

    def getNode(self, ID):
        return self.nodes.get(ID, None)

    def getEdgeBtwNodes(self, u, v):
        """
        :param u: a Node object
        :param v: a Node object
        :return: the Edge object connecting u and v if any, else None
        """
        r1 = repr(Edge(u.ID, v.ID))
        edge = self.edges.get(r1, None)
        if edge is None:
            r2 = repr(Edge(v.ID, u.ID))
            edge = self.edges.get(r2, None)
        return edge

    def getNumNodes(self):
        return len(self.nodes)

    def getNumEdges(self):
        return len(self.edges)

    def isEmpty(self):
        return self.getNumNodes() == 0

    def isTree(self):
        n = self.getNumNodes()
        m = self.getNumEdges()
        return m == n - 1

    def recomputeMaxDeg(self):
        self.maxDeg = max([N.degree for N in self.nodes.values()])

    def pushNodePoses(self):
        poses = {}
        for ID in self.nodes:
            node = self.nodes[ID]
            poses[ID] = (node.x, node.y)
        self.nodePosStack.append(poses)

    def popNodePoses(self):
        poses = self.nodePosStack.pop()
        for ID in self.nodes:
            node = self.nodes[ID]
            pos = poses[ID]
            node.x, node.y = pos

    def dropNodePoses(self):
        self.nodePosStack.pop()

    def fdlayout(self, iel, op, ccs=None, pcs=[], topoAddon=None, useNeighbourStress=False):
        alg, rs = self.setupFDLayout(iel, op=op, ccs=ccs, pcs=pcs, topoAddon=topoAddon,
                                     useNeighbourStress=useNeighbourStress)
        alg.run()
        self.moveNodesToRects(rs)

    def makeFeasible(self, op=False, ccs=None, pcs=[]):
        iel = 1 # doesn't matter
        alg, rs = self.setupFDLayout(iel, op=op, ccs=ccs, pcs=pcs, topoAddon=topoAddon)
        alg.makeFeasible()
        self.moveNodesToRects(rs)

    def setupFDLayout(self, iel, op=False, ccs=None, pcs=[], topoAddon=None, logger=None,
                      useNeighbourStress=False):
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('rs_es_ix')
        rs, es, ix = self.writeRsEsIx()
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('adg_build_FD_layout')
        alg = adg.ConstrainedFDLayout(rs, es, iel, op)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        alg.m_useNeighbourStress = useNeighbourStress

        if ccs is None:
            ccs = adg.CompoundConstraintPtrs()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('nodeconf_buildCCs')
        ncs = self.nodeConf.buildCCS(ix=ix)
        for nc in ncs:
            ccs.push_back(nc)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('extraPCs')
        for pc in self.extraPCs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        for pch in self.PCHolders:
            ccs2 = pch.getCCs(rs=rs, ix=ix)
            for cc in ccs2:
                ccs.push_back(cc)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        for pc in pcs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        alg.setConstraints(ccs)
        if topoAddon is not None:
            #topo = adg.ColaTopologyAddon()
            alg.setTopology(topoAddon)
            #alg.makeFeasible()
        return (alg, rs)

    def setupMajLayout(self, iel, op=False, ccs=None, pcs=[], logger=None,
                      useNeighbourStress=False, useScaling=False):
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('rs_es_ix')
        rs, es, ix = self.writeRsEsIx()
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('adg_build_Maj_layout')
        eLengths = adg.Doubles()
        for e in es:
            eLengths.push_back(1.0)
        alg = adg.simpleCMLFactory(rs, es, None, iel, useNeighbourStress)

        if op:
            alg.setAvoidOverlaps(True)
        alg.setScaling(useScaling)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        #alg.m_useNeighbourStress = useNeighbourStress

        if ccs is None:
            ccs = adg.CompoundConstraintPtrs()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('nodeconf_buildCCs')
        ncs = self.nodeConf.buildCCS(ix=ix)
        for nc in ncs:
            ccs.push_back(nc)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.startNewTimer('extraPCs')
        for pc in self.extraPCs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        for pch in self.PCHolders:
            ccs2 = pch.getCCs(rs=rs, ix=ix)
            for cc in ccs2:
                ccs.push_back(cc)
        if logger is not None and logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()

        for pc in pcs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        alg.setConstraintsVector(ccs)
        #alg.setConstraints(ccs)
        return (alg, rs)

    def listNodesAtDistanceFrom(self, n, node):
        rs, es, ix = self.writeRsEsIx()
        fdl = adg.ConstrainedFDLayout(rs, es, 1, False)
        d = fdl.getDistancesFromRect(ix(node))
        indices = []
        for i in range(len(d)):
            dist = d[i]
            if dist == n:
                indices.append(i)
        nodes = []
        for index in indices:
            for u in self.nodes.values():
                if ix(u) == index:
                    nodes.append(u)
        return nodes

    def shakeWithSolidEdges(self, iel, ccs=None, pcs=[], topoAddon=None, yFirst=False,
                            logger=None, gmlSuffix=None, useNeighbourStress=False, useScaling=False):
        """
        :param yFirst: say if you want to work in the y-dimension before the x-dimension
        Other parameters are as for fdlayout method.
        :return: nothing

        We do an fdlayout in each dimension, with overlap prevention and solidified
        axis-aligned edges.
        """
        axes = ['y', 'x'] if yFirst else ['x', 'y']
        for axis in axes:
            if logger is not None: logger.startT('copy_graph')
            H = self.copy()
            if logger is not None: logger.stopT()
            # For the moment want no extra gap, since this would be enforced even
            # between proper nodes and the edgenodes beside them, stretching things
            # out far more than necessary.
            if axis == 'x':
                H.nodeConf.extraGapY = 0
            else:
                H.nodeConf.extraGapX = 0
            # Now solidify the edges, and proceed.
            if logger is not None: logger.startT('solidify_edges')
            H.solidifyEdges(dirs=Compass.vertical if axis=='x' else Compass.horizontal)
            if logger is not None: logger.stopT()
            if logger is not None and gmlSuffix is not None:
                logger.writeGML(gmlSuffix+"_"+axis, graph=H)
            xAxis, yAxis = (True, False) if axis=='x' else (False, True)
            if useNeighbourStress:
                alg, rs = H.setupMajLayout(iel, op=True, logger=logger,
                                           useNeighbourStress=useNeighbourStress,
                                           useScaling=useScaling)
                #alg.m_doYAxisFirst = yAxis
                if logger is not None: logger.startT('adg_run_Maj_layout')
                alg.run(xAxis, yAxis)
                if logger is not None: logger.stopT()
                if logger is not None: logger.startT('move_nodes_to_rects')
                self.moveNodesToRects(rs)
                if logger is not None: logger.stopT()
            else:
                alg, rs = H.setupFDLayout(iel, op=True, logger=logger, useNeighbourStress=useNeighbourStress)
                alg.m_doYAxisFirst = yAxis
                # Set UCIs to see if anything went wrong.
                uciX = adg.UnsatisfiableConstraintInfoPtrs()
                uciY = adg.UnsatisfiableConstraintInfoPtrs()
                alg.setUnsatisfiableConstraintInfo(uciX, uciY)
                if logger is not None and logger.level >= LogLevel.DEBUG:
                    alg.outputInstanceToSVG("testOut/shake_%d_%s_init" % (axes.index(axis), axis))
                if logger is not None: logger.startT('adg_run_FD_layout')
                alg.run(xAxis, yAxis)
                if logger is not None: logger.stopT()
                if logger is not None and logger.level >= LogLevel.DEBUG:
                    alg.outputInstanceToSVG("testOut/shake_%d_%s_post" % (axes.index(axis), axis))
                if logger is not None: logger.startT('move_nodes_to_rects')
                self.moveNodesToRects(rs)
                if logger is not None: logger.stopT()
                if logger is not None and logger.level >= LogLevel.DEBUG:
                    print '\nShake x: %s, y: %s.\nUnsat-x: %d, Unsat-y: %d' % (
                        xAxis, yAxis, uciX.size(), uciY.size()
                    )

    def applyProjSeq(self, ps, logger, **kwargs):
        """
        Attempt to apply all the projections given by a ProjSeq object.
        Give up as soon as any of them fails.
        :param ps: a ProjSeq object
        :param kwargs: configuration options -- see code below
        :return: boolean saying whether all the projections were successful
        """
        # Set options:
        # A positive ideal edge length should be passed if stress is to be evaluated.
        iel = kwargs.get('iel', 0)
        # Overlap prevention
        opx = kwargs.get('opx', False)
        opy = kwargs.get('opy', False)
        # Edge solidification
        solidEdgesX = kwargs.get('solidEdgesX', False)
        solidEdgesY = kwargs.get('solidEdgesY', False)
        # Acceptance level:
        accept = kwargs.get('accept', 0)
        # Make accessible by dimension.
        opts = {
            adg.XDIM: {
                'op': opx, 'se': solidEdgesX
            },
            adg.YDIM: {
                'op': opy, 'se': solidEdgesY
            }
        }
        # Project
        if iel > 0:
            if logger.level >= LogLevel.TIMING:
                logger.startNewTimer('compute_stress')
            lastStress = self.computeStress(iel, logger=logger)
            if logger.level >= LogLevel.TIMING:
                logger.stopLastTimer()
        allOK = True
        for pcs, dim in ps:
            DEBUGGING = False
            if DEBUGGING:
                print 'Projection sequence set %s' % pcs
            if len(pcs) == 0: continue
            op, se = opts[dim]['op'], opts[dim]['se']
            result = self.project(logger, dim, op=op, pcs=pcs, solidEdges=se, accept=accept)
            if logger is not None:
                logger.addProj()
            if result > accept:
                allOK = False
                if DEBUGGING:
                    print '    Projection failed.'
            else:
                if DEBUGGING:
                    print '    Projection okay.'
                    logger.writeGML('lastSuccessfulProjection', graph=self)
                pass
            if iel > 0:
                if logger.level >= LogLevel.TIMING:
                    logger.startNewTimer('compute_stress')
                stress = self.computeStress(iel, logger=logger)
                if logger.level >= LogLevel.TIMING:
                    logger.stopLastTimer()
                dS = stress - lastStress
                ps.addStressChange(dS)
                lastStress = stress
            if not allOK:
                break
        return allOK

    def project(self, logger, dim, op=False, ccs=None, pcs=[], solidEdges=False,
                debugFilename=None, accept=0):
        """
        Project onto a set of constraints.
        :param dim: XDIM or YDIM -- the dimension in which to project
        :param op: boolean, saying whether you want overlap prevention constraints to be
                   generated
        :param ccs: any cola::CompoundConstraints you want to apply
        :param pcs: any PyConstraints you want to apply
        :param solidEdges: whether you want the axis-aligned edges to be solidified
                           before generating overlap prevention constraints
        :param debugFilename: (for debugging)
        :param accept: highest return value for which node positions will be accepted; see return
        :return: integer describing the result:
                    0: all constraints satisfiable,
                    1: some unsatisfiable, but only nonoverlap constraints
                    2: some unsatisfiable which are not nonoverlap constraints

        Node positions will be updated only if all constraints were satisfied.

        Note: we generate all constraints implied by the Graph's NodeConfig object.
        """
        if logger.level >= LogLevel.TIMING:
                logger.startNewTimer('set-up_projection')
        if logger.level >= LogLevel.DEBUG:
            if debugFilename is None: debugFilename = 'projTest'
            cpp = ''
        if solidEdges:
            H = self.copy()
            if logger.level >= LogLevel.TIMING:
                logger.startNewTimer('solidify_edges')
            H.solidifyEdges(dirs=Compass.vertical if dim==adg.XDIM else Compass.horizontal)
            if logger.level >= LogLevel.TIMING:
                logger.stopLastTimer()
        else:
            H = self
        if logger.level >= LogLevel.TIMING:
            logger.startNewTimer('rs_es_ix')
        rs, es, ix = H.writeRsEsIx()
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        if logger.level >= LogLevel.DEBUG:
            cpp += '/*\n'
            cpp += ' * Mapping of Node IDs to Variable IDs:\n'
            for ID in sorted(H.nodes.keys()):
                node = H.nodes[ID]
                cpp += ' * Node %d: Variable %d\n' % (ID, ix(node))
            cpp += '*/\n'
            with open('testOut/project/%s_solidEdges.gml' % debugFilename, 'w') as gmlFile:
                gmlFile.write(H.writeGML())
        n = rs.size()
        if logger.level >= LogLevel.DEBUG:
            cpp += H.writeCpp()
        if ccs is None:
            ccs = adg.CompoundConstraintPtrs()
        if logger.level >= LogLevel.TIMING:
            logger.startNewTimer('build_CCs')
        ncs = H.nodeConf.buildCCS(ix=ix)
        for nc in ncs:
            ccs.push_back(nc)
        for pc in self.extraPCs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        for pch in self.PCHolders:
            ccs2 = pch.getCCs(rs=rs, ix=ix)
            for cc in ccs2:
                ccs.push_back(cc)
        for pc in pcs:
            pcCCs = pc.buildCCs(rs=rs, ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
            if logger.level >= LogLevel.DEBUG:
                cpp += pc.writeCpp(ix=ix);
                cpp += 'ccs.push_back(%s);\n' % 'sep' if isinstance(pc, SepCo) else 'algn'
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        if logger.level >= LogLevel.DEBUG:
            cpp += 'cola::ProjectionResult result = cola::projectOntoCCs(%s, rs, ccs, %s, %d);\n' % (
                {adg.XDIM: 'vpsc::XDIM', adg.YDIM: 'vpsc::YDIM'}[dim],
                'true' if op else 'false',
                accept
            )
            cpp += 'std::cout << result.errorLevel << std::endl;\n'
            cpp += 'std::cout << result.unsatinfo << std::endl;\n'
            with open('testOut/project/%s.cpp' % debugFilename,'w') as cppFile:
                cppFile.write(cpp)
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
            logger.startNewTimer('adg_project')
        result = adg.projectOntoCCs(dim, rs, ccs, op, accept)
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        if result.errorLevel <= accept:
            self.moveNodesToRects(rs, ix=ix)
        elif logger.level >= LogLevel.DEBUG:
            print result.unsatinfo
            T = self.mapUnsatInfoToNodeLabels(H, result.unsatinfo, ix, dim)
            print T
            logger.writeGML('__UNSAT', graph=H, highlightNodes=[t[1] for t in T])
        return result.errorLevel

    def mapUnsatInfoToNodeLabels(self, H, unsatinfo, ix, dim,
                                 printLabelConstraints=True, printIDConstraints=False):
        """
        :param H: the graph
        :param unsatinfo: the string result.unsatinfo from result object returned by adg.projectOntoCCs
        :param ix: the usual node ID / rectangle index mapping
        :param dim: the dimension in which we were projecting
        :return: list of triples (rect_index, node_id, node_label)
        """
        # Invert the index map. We want a map from rect indices back to node IDs.
        ix2ID = {}
        for ID, node in H.nodes.items():
            i = ix(node.ID)
            ix2ID[i] = ID
        # Find all the named variables.
        lines = unsatinfo.split('\n')
        varnums = set()
        constraints = []

        class MiniConstraint:
            def __init__(self, vl, op, sep, rel, vr):
                self.vl=vl;self.op=op;self.sep=sep;self.rel=rel;self.vr=vr
            def findTriples(self, triples):
                I = [t[0] for t in triples]
                il, ir = I.index(int(self.vl)), I.index(int(self.vr))
                self.tl = triples[il]
                self.tr = triples[ir]
            def writeWithIDs(self, triples):
                self.findTriples(triples)
                return 'Node_%s %s %s %s Node_%s' % (
                    self.tl[1], self.op, self.sep, self.rel, self.tr[1]
                )
            def writeWithLabels(self, triples):
                self.findTriples(triples)
                return 'Node_%s %s %s %s Node_%s' % (
                    self.tl[2], self.op, self.sep, self.rel, self.tr[2]
                )

        for line in lines:
            if len(line) == 0 or line[0] != 'v': continue
            varnums.update(re.findall('v_(\d+)', line))
            try:
                vl, op, sep, rel, vr = re.match('v_(\d+) ([+-]) (\S+) (<=|==) v_(\d+)', line).groups()
                print '%s, %s' % (vl, vr)
                constraints.append(MiniConstraint(vl, op, sep, rel, vr))
            except:
                print 'Non-matching line:\n' + line

        varnums = list(sorted([int(n) for n in varnums]))

        ids = []
        for n in varnums:
            ID = ix2ID.get(n, None)
            if ID is None:
                ids.append(-1)
            else:
                ids.append(ID)

        labels = []
        for ID in ids:
            if ID == -1:
                labels.append('-1')
            else:
                labels.append(H.nodes[ID].label)

        triples = zip(varnums, ids, labels)
        if printLabelConstraints:
            for c in constraints:
                print c.writeWithLabels(triples)
        if printIDConstraints:
            for c in constraints:
                print c.writeWithIDs(triples)
        return triples

    def solidifyEdges(self, dirs=Compass.cwCards):
        """
        :param dirs: the Compass directions for which you want to generate nodes
        :return: nothing

        For each edge (u, v) for which u and v are configured in one of the direction
        in the dirs list, we add to G a node representing the edge, and configure it
        to stay between its endpoint nodes.
        """
        EDGETHICKNESS = Graph.EDGENODE_THICKNESS
        # Do not make rectangles that would be shorter than a certain epsilon length.
        # (If they are zero length, adaptagrams fails. And if they are nearly zero and
        # you write your graph as Dunnart SVG with only a few decimal places, they may
        # come out as literally zero.)
        MINLENGTH = 0.01
        GAP = 1
        for edge in self.edges.values():
            u, v = edge.src, edge.tgt
            d = self.nodeConf.getDirec(u, v)
            if d in dirs:
                node = Node()
                node.ID = self.getNextID()
                node.fill = "#4050C0"
                node.setIDAsLabel()
                ub, vb = u.boundingBoxxXyY(), v.boundingBoxxXyY()
                if d in Compass.vertical:
                    node.w = EDGETHICKNESS
                    node.h = min(
                        abs(vb[2]-ub[3]), abs(ub[2]-vb[3])
                    ) - 2*GAP
                    if node.h < MINLENGTH: continue
                    node.x = (u.x + v.x) / 2.0
                    node.y = min(ub[3], vb[3]) + node.h / 2.0 + GAP
                elif d in Compass.horizontal:
                    node.h = EDGETHICKNESS
                    node.w = min(
                        abs(vb[0]-ub[1]), abs(ub[0]-vb[1])
                    ) - 2*GAP
                    if node.w < MINLENGTH: continue
                    node.y = (u.y + v.y) / 2.0
                    node.x = min(ub[1], vb[1]) + node.w / 2.0 + GAP
                self.addNode(node)
                # NB: It is very important that we free the constraint on u and v in addition
                # to setting the new constraints on u and the new node, and on v and the new
                # node. Otherwise we will pass redundant equality constraints to VPSC, which
                # will report an unsatisfiability.
                self.nodeConf.free(u, v)
                self.nodeConf.setDirec(u, node, d)
                self.nodeConf.setDirec(node, v, d)

    def route(self, opts={}, iel=100):
        routing = opts.get('routing', adg.PolyLineRouting)
        router = adg.Router(routing)

        router.RubberBandRouting = opts.get('rubberBand', False)
        router.UseLeesAlgorithm = opts.get('Lees', True)
        router.InvisibilityGrph = opts.get('invis', True)

        # Hacky option I'm providing for certain test cases!
        fixRoutesBetweenRects = opts.get('fixRoutesBetweenRects', False)

        # Another slightly hacky option:
        svgOutputName = opts.get('svgOutputName', None)

        #router.setRoutingOption(adg.nudgeSharedPathsWithCommonEndPoint, False)
        router.setRoutingParameter(adg.crossingPenalty, opts.get('crossingPenalty', 0))
        #router.setRoutingParameter(adg.fixedSharedPathPenalty, 1000000)
        #router.setRoutingParameter(adg.segmentPenalty, 25)
        edgeReptoConnRef = {}
        maxNodeID = -1

        # We'll build a ColaTopologyAddon while we're at it.
        topoNodes = adg.TopologyNodePtrs()
        topoEdges = adg.TopologyEdgePtrs()
        rs, es, ix = self.writeRsEsIx()
        topoNodes.resize(rs.size())
        topoEdges.resize(es.size())

        for ID in self.nodes:
            node = self.nodes[ID]
            poly = node.libavoidPoly()
            sr = adg.ShapeRef(router, poly, ID)
            #nodeIDtoShapeRef[ID] = sr
            maxNodeID = max(maxNodeID, ID)
            index = ix(node)
            topoNodes[index] = adg.Node(index, rs[index])
        baseEdgeID = maxNodeID + 1
        for i, rep in enumerate(self.edges.keys()):
            edge = self.edges[rep]
            ID = baseEdgeID + i
            cr = adg.ConnRef(router, ID)
            edgeReptoConnRef[rep] = cr
            src, tgt = edge.libavoidConnEnds()
            cr.setEndpoints(src, tgt)
            #cr.setHateCrossings(True)
            if fixRoutesBetweenRects:
                ss, ts = edge.src.shape, edge.tgt.shape
                print [ss, ts]
                if ss == 'rectangle' and ts == 'rectangle':
                    print '    setting fixed route'
                    fr = adg.Polygon()
                    fr._id = ID
                    fr.ps.resize(2);
                    fr.ps[0] = adg.Point(edge.src.x, edge.src.y)
                    fr.ps[1] = adg.Point(edge.tgt.x, edge.tgt.y)
                    cr.setFixedRoute(fr)
        router.processTransaction()
        if svgOutputName is not None:
            router.outputInstanceToSVG('testOut/%s' % svgOutputName)
        # Store the connector routes in the Edge objects.
        for rep in edgeReptoConnRef:
            cr = edgeReptoConnRef[rep]
            # Put display route in Edge object
            polyline = cr.displayRoute()
            verbose = False
            if verbose:
                print 'Edge: %s' % rep
                print ' '*4 + ', '.join([
                    '(%.2f, %.2f)' % (p.x, p.y) for p in polyline.ps
                ])
            edge = self.edges[rep]
            edge.route = [
                [p.x, p.y] for p in polyline.ps
            ]
            # Make a topology route
            polyline = cr.route()
            eps = adg.TopologyEdgePointPtrs()
            eps.push_back(
                adg.EdgePoint(topoNodes[ix(edge.src)], adg.EdgePoint.CENTRE)
            )
            for j in range(1, polyline.size() - 1):
                p = polyline.ps[j]
                eps.push_back(
                    adg.EdgePoint(
                        topoNodes[ix(self.nodes[p.id])],
                        {
                            0: adg.EdgePoint.BR,
                            1: adg.EdgePoint.TR,
                            2: adg.EdgePoint.TL,
                            3: adg.EdgePoint.BL
                        }.get(p.vn, adg.EdgePoint.CENTRE)
                    )
                )
            eps.push_back(
                adg.EdgePoint(topoNodes[ix(edge.tgt)], adg.EdgePoint.CENTRE)
            )
            i = topoEdges.size()
            topoEdges.push_back(
                adg.Edge(i, iel, eps)
            )
        # Build and return topology addon.
        topo = adg.ColaTopologyAddon(topoNodes, topoEdges)
        return topo

    def inferNodeConf(self):
        """
        Set self.nodeConf to a fresh NodeConfig object enforcing precisely those
        configurations that are already present in the graph.
        :return: nothing
        """
        self.nodeConf = NodeConfig(self)
        for edge in self.edges.values():
            n1, n2 = edge.src, edge.tgt
            if n1.x == n2.x or n1.y == n2.y:
                direc = Compass.cardinalDirection(n1, n2)
                #print 'Setting direction %d from node %2d to node %2d' % (
                #    direc, n1.ID, n2.ID
                #)
                self.nodeConf.setDirec(n1, n2, direc)

    def writeRsEsIx(self):
        """
        Returns triple (rs, es, ix) being adaptagrams rectanges and edges
        to represent this graph, and a map (callable) from Node objects
        to their indices in the rs vector.
        """
        rs = adg.RectanglePtrs()
        es = adg.ColaEdges()
        L = sorted(self.nodes.keys())
        ix = lambda N: L.index(N)
        for ID in L:
            N = self.nodes[ID]
            x,X,y,Y = N.getxXyY()
            R = adg.Rectangle(x,X,y,Y)
            rs.push_back(R)
        for rep in self.edges:
            E = self.edges[rep]
            six = ix(E.src.ID); tix = ix(E.tgt.ID)
            E = adg.ColaEdge(six, tix)
            es.push_back(E)
        return (rs, es, ix)

    def writeCpp(self):
        """
        :return: c++ code to build this graph
        """
        rs, es, ix = self.writeRsEsIx()
        s = ''
        s += 'vpsc::Rectangles rs;\n'
        s += 'vpsc::Rectangle *rect;\n'
        s += 'std::vector<cola::Edge> es;\n'
        s += 'cola::Edge edge;\n'
        for r in rs:
            s += 'rect = new vpsc::Rectangle(%.2f, %.2f, %.2f, %.2f);\n'%(
                r.getMinX(), r.getMaxX(), r.getMinY(), r.getMaxY()
            )
            s += 'rs.push_back(rect);\n'
        for e in es:
            s += 'edge = cola::Edge(%d, %d);\n'%(
                e[0], e[1]
            )
            s += 'es.push_back(edge);\n'
        s += self.nodeConf.writeCpp(ix=ix)
        return s


    def moveNodesToRects(self, rs, ix=None):
        """
        Set node positions based on rectangle positions.
        For use after running a CoLa layout.
        """
        if ix is None:
            for R, ID in zip(rs, sorted(self.nodes.keys())):
                N = self.nodes[ID]
                N.moveToRect(R)
        else:
            for node in self.nodes.values():
                node.moveToRect(rs[ix(node.ID)])

    def computeAvgNodeDim(self):
        """
        Compute the average of all node heights and widths.
        """
        s = 0
        n = 0
        for N in self.nodes.values():
            s += N.w + N.h
            n += 2
        avg = s/float(n)
        return avg

    def getChainsAndCycles(self):
        """
        Identify all sequences of consecutive "links" (degree-2 nodes) in this graph.
        :return: dictionary of the form:
            {
                'chains': [
                    [nodes], [nodes], ..., [nodes]
                ],
                'cycles': [
                    [nodes], [nodes], ..., [nodes]
                ]
            }
        """
        chains = []
        cycles = []
        # Build /list/ of all links in the graph.
        # Used to use a set here instead of a list, but it is generally agreed upon:
        #     http://stackoverflow.com/questions/3406341/iteration-order-of-sets-in-python
        # that counting on a consistent iteration order over a set is a bad idea.
        # We sacrifice a modicum of speed in an effort to eradicate all randomness.
        allLinks = []
        for ID in sorted(self.nodes.keys()):
            node = self.nodes[ID]
            if node.degree == 2:
                allLinks.append(node)
        while len(allLinks) > 0:
            # Take next link from the list.
            # We take from the back of the list since we have the pop method for that.
            L0 = allLinks.pop()
            # Initialise deque with L0, to hold all links in
            # the chain to which L0 belongs.
            links = Deque([L0])
            # Get the two edges of L0, and prepare to explore in both directions.
            E0 = L0.edges.values()
            assert(len(E0)==2)
            direc = -1
            polygon = False
            for e in E0:
                if polygon: break
                # Explore from link L0 in one direction.
                direc = -direc
                last = L0
                done = False
                while not done:
                    # Consider the next node in the current direction.
                    next = e.otherEnd(last)
                    if next == L0:
                        # In this case the entire connected component to which L0
                        # belongs is a mere polygon.
                        polygon = True
                        cycles.append(links)
                        links = []
                        done = True
                    elif next.degree == 2:
                        # This must be a link which we have not encountered before.
                        allLinks.remove(next)
                        if direc == 1:
                            links.append(next)
                        elif direc == -1:
                            links.appendleft(next)
                        e = filter(lambda e1: e1 is not e, next.edges.values())[0]
                        last = next
                    else:
                        # We've reached the "anchor node" at one end of the chain.
                        done = True
            # Now have explored from link L0 in both directions, or else found that
            # it belonged to a polygon.
            if len(links) > 0:
                chains.append(links)
        return {'chains': chains, 'cycles': cycles}

    def connComps(self):
        """
        Return a list of Graphs containing all the connected
        components of this graph.

        The nodes and edges in the component graphs are /the same objects/
        as the nodes and edges in this graph.
        """
        C = []
        # Make a copy of the node dict to manage the process, but its
        # keys will point to the SAME node objects, NOT copies.
        nodes = {}
        for ID in self.nodes: nodes[ID] = self.nodes[ID]
        while nodes:
            H = Graph()
            # Start with any remaining node.
            k0 = nodes.keys()[0]
            u = nodes[k0]
            del nodes[k0]
            H.nodes[k0] = u
            # We will use a queue of pairs (e,u), being an edge e
            # along with the endpoint u that nominated it.
            el = u.edgeList()
            queue = [(e,u) for e in el]
            # BFS:
            while queue:
                e,u = queue[0]
                del queue[0]
                v = e.otherEnd(u)
                if H.nodes.has_key(v.ID): continue
                H.nodes[v.ID] = v
                H.edges[repr(e)] = e
                del nodes[v.ID]
                el = v.edgeList()
                queue.extend([(f,v) for f in el])
            H.recomputeMaxDeg()
            C.append(H)
        return C

    def severEdge(self,edge):
        edge.sever()
        del self.edges[repr(edge)]

    def deleteNode(self,node):
        del self.nodes[node.ID]

    def severNodes(self,nodes,buckets=None):
        """
        nodes: {ID:Node} dict of nodes to be removed from the
               graph. 
        buckets: optional NodeBuckets object which should be informed of
                 any changing node degrees.
        """
        for ID in nodes:
            u = nodes[ID]
            el = u.edgeList()
            for e in el:
                v = e.otherEnd(u)
                deg = v.degree
                self.severEdge(e)
                if buckets:
                    buckets.moveNode(v.ID,deg,deg-1)
                    # We don't move u in the buckets, since we
                    # presume that it has already been removed.
            self.nodeConf.removeNode(u)
            del self.nodes[u.ID]

    def identifyRootNode(self):
        """
        For use with trees which have been built using the treeSerialNo
        attribute of the nodes involved.
        """
        M = 0
        R = None
        for ID in self.nodes:
            u = self.nodes[ID]
            #sys.stderr.write('serial: %d\n'%u.serial)
            if u.treeSerialNo > M:
                M = u.treeSerialNo
                R = u
        self.rootNode = R
        R.setIsRootNode(True)

    def writeGML(self, idt='', highlightNodes=None):
        s = ''
        s += idt + 'graph\n[\n'
        s += self.writeNodesAndEdgesGML(idt=idt, highlightNodes=highlightNodes)
        s += idt + ']\n'
        return s

    def writeJSON(self):
        d = {}
        nodes = {}
        for node in self.nodes.values():
            nodes[node.ID] = node.buildJSONDict()
        d['nodes'] = nodes
        edges = {}
        for dsc, edge in self.edges.items():
            edges[repr(edge)] = edge.buildJSONDict()  # TODO
        d['edges'] = edges
        j = json.dumps(d, indent=4)
        return j

    def writeNodesAndEdgesGML(self, idt='', highlightNodes=None):
        s = ''
        for ID in sorted(self.nodes.keys()):
            N = self.nodes[ID]
            s += N.writeGML(idt=idt+'    ', highlightNodes=highlightNodes)
        for rep in sorted(self.edges.keys()):
            E = self.edges[rep]
            s += E.writeGML(idt=idt+'    ')
        return s

    def writeSVG(self, w, h, pad=40, randomise=False, rednodeIDs=[], ownFill=False):
        if randomise:
            nums = range(len(self.nodes))
            random.shuffle(nums)
            for i, node in enumerate(self.nodes.values()):
                node.label = str(nums[i])
        x0, y0, w0, h0 = self.boundingBoxXYWH(includeBends=True)
        FIT = False
        ZOOM = True
        scale = 1.0
        if FIT:
            w, h = w0 + pad, h0 + pad
        elif ZOOM:
            scale = min((w-pad)/w0, (h-pad)/h0)
        s = ''
        s += '<?xml version="1.0" encoding="UTF-8"?>\n'
        s += '<svg id="canvasSVG" xmlns="http://www.w3.org/2000/svg" '
        s += 'width="%spx" height="%spx">\n' % (w, h)
        #s += 'width="100%" height="100%">\n'
        # Background rectangle:
        s += svg.rect(0, 0, w, h, {'fill': '#bdbddf', 'stroke': 'none'}) + '\n'
        if ZOOM:
            w0, h0 = scale*w0, scale*h0
        dx, dy = (w-w0)/2.0, (h-h0)/2.0
        s += '<g id="outerGroup">\n'
        s += '<g transform="translate(%.2f,%.2f),scale(%.2f),translate(%.2f,%.2f)">\n' % (
            dx, dy, scale, -x0, -y0
        )
        for edge in self.edges.values():
            s += ' '*4 + edge.writeSVG() + '\n'
        for node in self.nodes.values():
            red = False
            if node.ID in rednodeIDs:
                red=True
            s += ' '*4 + node.writeSVG(red=red, ownFill=ownFill) + '\n'
        s += '</g>\n'
        s += '</g>\n'
        s += '</svg>\n'
        return s

    def writeDunnartSVG(self, opts=None):
        """
        :param opts: optional dictionary of Dunnart options, e.g.
                     { automaticGraphLayout:1 }
        """
        s = ''
        s += '<?xml version="1.0" encoding="UTF-8"?>\n'
        s += '<svg xmlns:dunnart="http://www.dunnart.org/ns/dunnart">\n'
        if opts is not None:
            s += '<dunnart:options '
            s += ' '.join([
                '%s="%s"' % (k, opts[k]) for k in opts
            ])
            s += '/>\n'
        maxID = -1
        for node in self.nodes.values():
            s += ' '*4 + node.writeDunnartSVG() + '\n'
            maxID = max(maxID, node.ID)
        nextID = maxID + 1
        for edge in self.edges.values():
            s += ' '*4 + edge.writeDunnartSVG(nextID) + '\n'
            nextID += 1
        #s += self.pcs.writeDunnartSVG(nextID, self)
        s += self.nodeConf.writeDunnartSVG(nextID)
        s += '</svg>'
        return s


class Components:
    """
    For putting several graphs together into one diagram.
    """

    def __init__(self, comps):
        self.comps = comps

    def addComps(self, comps):
        self.comps.extend(comps)

    def writeGML(self, idt=''):
        s = ''
        s += idt + 'graph\n[\n'
        #s += idt + '    ' + 'directed    1\n'
        maxID = 0
        for G in self.comps:
            minID, Gmax = G.getExtremeIDs()
            shift = max(0, maxID - minID + 1)
            G.shiftIDs(shift)
            maxID = Gmax + shift
            s += G.writeNodesAndEdgesGML(idt=idt)
        s += idt + ']\n'
        return s

######################################################################

def gml2dunnartShapes(gmlShapeName):
    return {
        'rectangle': 'rect',
        'circle': 'ellipse'
    }[gmlShapeName]

def newBendNode(ID, x, y):
    node = Node()
    node.ID = ID
    node.x, node.y = x, y
    node.w, node.h = 8, 8
    node.shape = 'circle'
    node.fill = '#00FFFF'
    return node


class RoutingRig:

    def __init__(self, opts):
        """
        :param opts: dictionary containing router options
        """
        routing = opts.get('routing', adg.PolyLineRouting)
        router = adg.Router(routing)
        self.router = router
        # Set options
        self.setOpts(opts)
        # Nodes and Edges
        self.nodes = {}
        self.edges = {}
        # Maps
        self.nodeIDtoShapeRef = {}
        self.edgeReptoConnRef = {}
        # For internal IDs for the ShapeRefs and ConnRefs:
        self.nextID = 0

    def setOpts(self, opts, router=None):
        if router is None:
            router = self.router
        parameters = [
            'segmentPenalty',
            'anglePenalty',
            'crossingPenalty',
            'clusterCrossingPenalty',
            'fixedSharedPathPenalty',
            'portDirectionPenalty',
            'shapeBufferDistance',
            'idealNudgingDistance',
            'reverseDirectionPenalty'
        ]
        options = [
            'nudgeOrthogonalSegmentsConnectedToShapes',
            'improveHyperedgeRoutesMovingJunctions',
            'penaliseOrthogonalSharedPathsAtConnEnds',
            'nudgeOrthogonalTouchingColinearSegments',
            'performUnifyingNudgingPreprocessingStep',
            'improveHyperedgeRoutesMovingAddingAndDeletingJunctions',
            'nudgeSharedPathsWithCommonEndPoint'
        ]
        fields = [
            'RubberBandRouting',
            'UseLeesAlgorithm',
            'InvisibilityGrph'
        ]
        noAction = [
            'routing'
        ]
        for k, v in opts.items():
            if k in parameters:
                a = adg.__dict__[k]
                router.setRoutingParameter(a, v)
            elif k in options:
                a = adg.__dict__[k]
                router.setRoutingOption(a, v)
            elif k in fields:
                router.__setattr__(k, v)
            elif k in noAction:
                pass
            else:
                raise Exception("RoutingRig did not recognise option: %s" % k)

    def takeNextID(self):
        "For internal use only."
        n = self.nextID
        self.nextID = n + 1
        return n

    def route(self, setRoutesInEdges=False):
        """
        Route the connectors.
        :param setRoutesInEdges: Set True if you want the connector routes to be stored
                                 in the Edge objects' self.route fields. Otherwise it's
                                 up to you to use the routes manually.
        :return: nothing
        """
        self.router.processTransaction()
        if setRoutesInEdges:
            for rep in self.edgeReptoConnRef:
                cr = self.edgeReptoConnRef[rep]
                polyline = cr.displayRoute()
                edge = self.edges[rep]
                edge.route = [
                    [p.x, p.y] for p in polyline.ps
                ]

    def getPtsFromConnRef(self, cr):
        """
        :param cr: a ConnRef object
        :return: list of points given as ordered pairs: [(x0, y0), ..., (xn, yn)]
        """
        return [(p.x, p.y) for p in cr.displayRoute().ps]

    def addConnectedPts(self, p0, p1, cd0=adg.ConnDirAll, cd1=adg.ConnDirAll):
        """
        Add a pair of connected points to the router.
        :param p0: point (x0, y0)
        :param p1: point (x1, y1)
        :return: the ConnRef object for the connector
        """
        pt0 = adg.Point(p0[0], p0[1])
        pt1 = adg.Point(p1[0], p1[1])
        end0 = adg.ConnEnd(pt0, cd0)
        end1 = adg.ConnEnd(pt1, cd1)
        crID = self.takeNextID()
        cr = adg.ConnRef(self.router, crID)
        cr.setEndpoints(end0, end1)
        return cr

    def addPolygonByPointList(self, pts):
        """
        :param pts: list of points (x, y) of a polygon to be added as obstacle
        :return: the ShapeRef object that is constructed and added
        """
        n = len(pts)
        poly = adg.Polygon(n)
        for i, pt in enumerate(pts):
            x, y = pt
            poly.setPoint(i, adg.Point(x, y))
        srID = self.takeNextID()
        sr = adg.ShapeRef(self.router, poly, srID)
        return sr

    def addNodesAndEdges(self, nodes, edges, fixedStraightRoutes=False, connDirs={}):
        """
        Add nodes and edges to the router.
        :param nodes: dictionary of the form ID:Node (may be empty)
        :param edges: dictionary of the form dsc:Edge (may be empty)
        :param fixedStraightRoutes: Set True if you want each edge to be given a fixed
                                    straight route from endpt to endpt. Presumably you
                                    will add other edges that don't have fixed routes,
                                    so the router has something to do!
        :param connDirs: If desired, pass a dictionary whose keys equal those of the
                         edges dictionary, and where each value is a pair (ds, dt)
                         giving the libavoid ConnDir flags for the src and tgt ends
                         of the corresponding edge, respectively.
        :return: nothing
        """
        n2s = self.nodeIDtoShapeRef
        e2c = self.edgeReptoConnRef
        router = self.router
        for ID, node in nodes.items():
            self.nodes[ID] = node
            poly = node.libavoidPoly()
            srID = self.takeNextID()
            sr = adg.ShapeRef(router, poly, srID)
            n2s[ID] = sr
        for rep, edge in edges.items():
            self.edges[rep] = edge
            crID = self.takeNextID()
            cr = adg.ConnRef(router, crID)
            ds, dt = connDirs.get(rep, (adg.ConnDirAll, adg.ConnDirAll))
            src, tgt = edge.libavoidConnEnds(srcDirs=ds, tgtDirs=dt)
            cr.setEndpoints(src, tgt)
            e2c[rep] = cr
            if fixedStraightRoutes:
                fr = adg.Polygon()
                fr._id = crID
                fr.ps.resize(2);
                fr.ps[0] = adg.Point(edge.src.x, edge.src.y)
                fr.ps[1] = adg.Point(edge.tgt.x, edge.tgt.y)
                cr.setFixedRoute(fr)

    def addGraphToRouter(self, G, fixedStraightRoutes=False):
        """
        Add all the nodes and edges of a Graph to self.router.
        :param G: a Graph object
        :return: nothing

        See addNodesAndEdges method.
        """
        nodes = G.nodes
        edges = G.edges
        self.addNodesAndEdges(nodes, edges, fixedStraightRoutes=fixedStraightRoutes)


class Node:

    def __init__(self, intcoords=False):
        self.__dict__['intcoords'] = intcoords
        self.graph = None
        self.ID = 0
        self.edges = {}
        self.degree = 0
        self.treeSerialNo = 0
        self.label = ""
        if intcoords:
            self.x, self.y, self.w, self.h = 0, 0, 30, 30
        else:
            self.x, self.y, self.w, self.h = 0.0, 0.0, 30.0, 30.0
        self.shape = 'rectangle'
        self.fill = '#CCCCEE' #(Dunnart node colour)
        self.stroke = '#000'
        self.isRootNode = False
        self.useYEdColour = False
        self.linkFill = '#008000' # special colour for degree-2 nodes
        self.rootFill = '#F08010' # special colour for root nodes
        self.yEdFill = '#F0C000'
        self.highlightFill = '#F00000'
        self.lowlightFill = '#F0C000'

    def __repr__(self):
        s = ''
        s += 'Node %s: (x,y)=(%.2f,%.2f) (w,h)=(%.1f,%.1f)'%(
            self.ID, self.x, self.y, self.w, self.h
        )
        if self.isRootNode:
            s += ' *'
        return s

    def addPadding(self, xPad, yPad):
        self.w += xPad
        self.h += yPad

    def buildJSONDict(self):
        d = {
            'ID': self.ID,
            'x': self.x, 'y': self.y,
            'w': self.w, 'h': self.h
        }
        return d

    def __setattr__(self, key, value):
        d = self.__dict__
        if self.intcoords and key in ['x', 'y', 'w', 'h']:
            d[key] = int(round(value))
        else:
            d[key] = value

    def centre(self):
        return (self.x, self.y)

    def liesOppositeSegment(self, seg, openInterval=False):
        """
        :param seg: a LineSegment object
        :param openInterval: set True if you want to know whether this node
                             lies opposite the segment's /open/ interval, not closed
        :return: boolean saying whether this Node's bbox lies beside the segment
        """
        DEBUG = False
        x, X, y, Y = self.boundingBoxxXyY()
        I = [x, X] if seg.varDim == adg.XDIM else [y, Y]
        if openInterval:
            ans = seg.openIntervalIntersects(I)
        else:
            ans = seg.closedIntervalIntersects(I)
        if DEBUG:
            print '    Node %s with %s-interval [%s, %s] %s lie opposite segment\n        %s' % (
                self.ID, {adg.XDIM: 'x', adg.YDIM: 'y'}[seg.varDim],
                I[0], I[1],
                'DOES' if ans else 'does NOT',
                seg
            )
            if ans:
                print '    Therefore:'
                print '        %s > %s: %s' % (
                    I[1], seg.wl, (I[1] > seg.wl)
                )
                print '    and:'
                print '        %s > %s: %s' % (
                    seg.wh, I[0], (seg.wh > I[0])
                )
        return ans

    def somePointOppositeSegment(self, seg):
        x, X, y, Y = self.boundingBoxxXyY()
        I = [x, X] if seg.varDim == adg.XDIM else [y, Y]
        J = seg.closedIntervalIntersection(I)
        if J is None:
            return None
        else:
            w = J[0]
            pt = (w, self.y) if seg.varDim == adg.XDIM else (self.x, w)
            return pt

    def closedIntervalOppositeSegment(self, seg):
        x, X, y, Y = self.boundingBoxxXyY()
        I = [x, X] if seg.varDim == adg.XDIM else [y, Y]
        J = seg.closedIntervalIntersection(I)
        return J

    def halfDims(self):
        if self.intcoords:
            hw = int(ceil(self.w/2.0))
            hh = int(ceil(self.h/2.0))
        else:
            hw, hh = self.w/2.0, self.h/2.0
        return (hw, hh)

    def getBdryCompassPt(self, direc):
        """
        :param direc: a Compass direction
        :return: point (u, v) on the boundary of this node in the given direction
                 from the centre
        """
        cx, cy = self.x, self.y
        hw, hh = self.halfDims()
        sgnx, sgny = Compass.vectorSigns(direc)
        u, v = cx + sgnx * hw, cy + sgny * hh
        return (u, v)

    def boundingBoxxXyY(self):
        """
        :return: bounding box in the form (x, X, y, Y) giving extreme coords
        """
        hw, hh = self.halfDims()
        u, v = self.x - hw, self.y - hh
        # Now use 2*hw, 2*hh instead of self.w, self.h, since in integer
        # case the halfdims are rounded up with the ceiling function.
        return (u, u + 2*hw, v, v+ 2*hh)

    def boundingBoxXYWH(self):
        """
        :return: bounding box of the node as (ULCx, ULCy, W, H)
        """
        hw, hh = self.halfDims()
        u, v = self.x - hw, self.y - hh
        return (u, v, 2*hw, 2*hh)

    def setIDAsLabel(self):
        self.label = str(self.ID)

    def isLink(self):
        return self.degree == 2

    def getFill(self, highlightNodes=None):
        if highlightNodes is not None:
            return self.highlightFill if self.ID in highlightNodes else self.lowlightFill
        elif self.useYEdColour:
            return self.yEdFill
        elif self.isRootNode and self.rootFill is not None:
            return self.rootFill
        elif self.isLink() and self.linkFill is not None:
            return self.linkFill
        else:
            return self.fill

    def writeSVG(self, red=False, ownFill=False):
        s = ''
        fill = '#f00' if red else '#f0f0d2'
        if ownFill:
            fill = self.fill
        s += svg.rectC(self.x, self.y, self.w, self.h, {
            'fill': fill,
            'onclick': 'ow3.nodeClick(this)',
            'id': str(self.ID)
        }) + '\n'
        #s += svg.text(self.x, self.y, self.label, {'centre': True}) + '\n'
        return s

    def writeDunnartSVG(self):
        s = ''
        s += '<dunnart:node type="%s" id="%d" ' % (
            gml2dunnartShapes(self.shape), self.ID
        )
        s += 'cx="%.2f" cy="%.2f" ' % (self.x, self.y)
        s += 'width="%.2f" height="%.2f" ' % (self.w, self.h)
        s += 'fillColour="%sFF" label="%s"/>' % (
            self.getFill()[1:], self.label
        )
        return s

    def libavoidPoly(self):
        poly = adg.Polygon(4)
        ulx = self.x - self.w/2.0
        uly = self.y - self.h/2.0
        poly.setPoint(0, adg.Point(ulx, uly))
        poly.setPoint(1, adg.Point(ulx + self.w, uly))
        poly.setPoint(2, adg.Point(ulx + self.w, uly + self.h))
        poly.setPoint(3, adg.Point(ulx, uly + self.h))
        return poly

    def getChildren(self):
        """
        :return: a list of the neighbours of this node that sit as
                 the target end of the connecting edge
        """
        kids = []
        for E in self.edges.values():
            if self is E.src:
                kids.append(E.tgt)
        return kids

    def getNeighbours(self):
        """
        :return: a list of all neighbours of this node
        """
        return [
            e.otherEnd(self) for e in self.edges.values()
        ]

    def getNbrsCWCyclic(self):
        """
        :return: list of all neighbours in clockwise cyclic order
                (assuming the usual graphics convention of x increasing
                 to the right and y increasing downward)
        """
        nbrs = self.getNeighbours()
        nbrs.sort(key=lambda n: math.atan2(
            n.y - self.y, n.x - self.x
        ))
        return nbrs

    def moveToRect(self, R):
        self.x = R.getCentreX()
        self.y = R.getCentreY()

    def setIsRootNode(self, b):
        self.isRootNode = b
        if b:
            self.fill = '#FF0000'

    def excisedCopy(self):
        """
        Return a shallow copy of this node, with graph, edges, and degree
        reset, as though this node had been excised from its graph.
        """
        C = copy.copy(self)
        C.graph = None
        C.edges = {}
        C.degree = 0
        return C

    def setGraph(self, G):
        self.graph = G

    def initFromGML(self, d):
        self.ID = d['id']
        self.label = d.get('label','')
        g = d['graphics']
        self.x = g['x']
        self.y = g['y']
        self.w = g['w']
        self.h = g['h']
        self.shape = g['type']
        self.fill = g.get('fill', '#FFCC00')
        self.stroke = g.get('outline', '#000000')

    def addEdge(self, E):
        self.edges[repr(E)] = E
        self.degree += 1

    def loseEdge(self, E):
        del self.edges[repr(E)]
        self.degree -= 1

    def edgeList(self):
        return self.edges.values()

    def getxXyY(self):
        x = self.x - self.w/2.0
        X = x + self.w
        y = self.y - self.h/2.0
        Y = y + self.h
        return (x,X,y,Y)

    def writeGML(self, idt='', highlightNodes=None):
        s = ''
        I = idt + '    '
        s += idt + 'node\n'
        s += idt + '[\n'
        s += I + 'id    %d\n'%self.ID
        s += I + 'label    "%s"\n'%self.label
        # Graphics
        s += I + 'graphics\n'
        s += I + '[\n'
        s += I + '    x    %.2f\n'%self.x
        s += I + '    y    %.2f\n'%self.y
        s += I + '    w    %.2f\n'%self.w
        s += I + '    h    %.2f\n'%self.h
        s += I + '    type    "%s"\n'%self.shape
        s += I + '    fill    "%s"\n'%self.getFill(highlightNodes=highlightNodes)
        s += I + '    outline    "%s"\n'%self.stroke
        s += I + ']\n'
        # End graphics
        s += idt + ']\n'
        return s



######################################################################

class Edge:

    def __init__(self, srcID=0, tgtID=0):
        self.srcID = srcID
        self.tgtID = tgtID
        self.src = None
        self.tgt = None
        self.fill = '#000000'
        self.route = None
        self.routePoints = []
        self.bends = []

    def buildJSONDict(self):
        d = {
            'source': self.src.ID, 'target': self.tgt.ID,
            'route': self.route
        }
        return d

    def translateBy(self, vect):
        if self.route is None: return
        dx, dy = vect
        route = []
        for p in self.route:
            x, y = p
            q = (x + dx, y + dy)
            route.append(q)
        self.route = route

    def rotate90cw(self):
        if self.route is None: return
        route = []
        for p in self.route:
            x, y = p
            q = (-y, x)
            route.append(q)
        self.route = route

    def rotate90acw(self):
        if self.route is None: return
        route = []
        for p in self.route:
            x, y = p
            q = (y, -x)
            route.append(q)
        self.route = route

    def rotate180(self):
        if self.route is None: return
        route = []
        for p in self.route:
            x, y = p
            q = (-x, -y)
            route.append(q)
        self.route = route

    def clearRoute(self):
        self.route = None
        self.routePoints = []

    def addBendPt(self, p):
        self.bends.append(p)

    def buildRoute(self):
        # Only build if have bends.
        if len(self.bends) == 0: return
        route = []
        route.append(self.src.centre())
        route.extend(self.bends)
        route.append(self.tgt.centre())
        self.route = route

    def libavoidConnEnds(self, srcDirs=adg.ConnDirAll, tgtDirs=adg.ConnDirAll):
        srcPt = adg.Point(self.src.x, self.src.y)
        tgtPt = adg.Point(self.tgt.x, self.tgt.y)
        srcEnd = adg.ConnEnd(srcPt, srcDirs)
        tgtEnd = adg.ConnEnd(tgtPt, tgtDirs)
        return (srcEnd, tgtEnd)

    def initFromGML(self, d):
        self.srcID = d['source']
        self.tgtID = d['target']
        # Check for route points.
        g = d['graphics']
        if g.has_key('Line'):
            pts = g['Line']['Line']
            pts = [[p['x'], p['y']] for p in pts]
            self.route = pts

    def setGraph(self, G):
        self.graph = G
        # Assume that by this time the Edge knows its srcID and tgtID,
        # and the graph has those nodes.
        self.src = G.getNode(self.srcID)
        self.tgt = G.getNode(self.tgtID)
        # Let the nodes know about the edge, and adjust their degree
        # accordingly.
        self.src.addEdge(self)
        self.tgt.addEdge(self)

    def __repr__(self):
        s = self.src.ID if self.src is not None else self.srcID
        t = self.tgt.ID if self.tgt is not None else self.tgtID
        return '%s --> %s'%(s, t)

    def otherEnd(self,node):
        return self.tgt if node.ID==self.src.ID else self.src

    def sever(self):
        self.src.loseEdge(self)
        self.tgt.loseEdge(self)

    def boundingBoxxXyY(self):
        #route = self.route[:]
        sc = self.src.centre()
        tc = self.tgt.centre()
        x = min(sc[0], tc[0])
        X = max(sc[0], tc[0])
        y = min(sc[1], tc[1])
        Y = max(sc[1], tc[1])
        if self.route is not None:
            for pt in self.route:
                u, v = pt
                x, X = min(u, x), max(u, X)
                y, Y = min(v, y), max(v, Y)
        return (x, X, y, Y)

    def writeGML(self, idt=''):
        s = ''
        I = idt + '    '
        s += idt + 'edge\n'
        s += idt + '[\n'
        s += I + 'source    %d\n'%self.srcID
        s += I + 'target    %d\n'%self.tgtID
        # Graphics
        s += I + 'graphics\n[\n'
        s += I + '    smoothBends    1\n'
        s += I + '    fill    "%s"\n'%self.fill
        if self.route:
            # GML is finicky about routes. The first and last points must equal
            # the centres of the endpt nodes. (Or maybe it is just yEd?)
            # So make sure this is the case. (libavoid does not ensure this.)
            route = self.route[:]
            sc = self.src.centre()
            tc = self.tgt.centre()
            if route[0] != sc:
                route.insert(0, sc)
            if route[-1] != tc:
                route.append(tc)
            s += '    Line\n    [\n'
            for pt in route:
                s += ' '*8 + 'point\n' + ' '*8 + '[\n'
                s += ' '*12 + 'x    %.2f\n' % pt[0]
                s += ' '*12 + 'y    %.2f\n' % pt[1]
                s += ' '*8 + ']\n'
            s += '    ]\n'
        s += I + ']\n'
        # End graphics
        s += idt + ']\n'
        return s

    def writeSVG(self):
        pts = self.route
        if pts is None:
            pts = [self.src.centre(), self.tgt.centre()]
        if len(pts) > 2:
            # A little fix for yEd-generated GML:
            # (Sometimes the endpt is not axis-aligned with the next point.)
            jobs = [(0, 1), (-1, -2)]
            for job in jobs:
                endpt = pts[job[0]]
                refpt = pts[job[1]]
                ex, ey = endpt
                rx, ry = refpt
                if ex != rx and ey != ry:
                    dx, dy = abs(ex - rx), abs(ey - ry)
                    if dx > dy:
                        endpt = [ex, ry]
                    else:
                        endpt = [rx, ey]
                    pts[job[0]] = endpt
        #print 'src: %d, tgt: %d, pts: %s' % (self.srcID, self.tgtID, pts)
        d = svg.roundedOrthoConnectorData(pts, curveRadius=10)
        #print '        %s' % d
        s = svg.path(d)
        return s

    def writeDunnartSVG(self, ID):
        s = ''
        s += '<dunnart:node type="connector" id="%d" srcID="%d" dstID="%d"/>'%(
            ID, self.srcID, self.tgtID
        )
        return s

