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

import adaptagrams.adaptagrams as adg
from gmlparse import buildGraph
from treePruner3 import prune
from graphs import RoutingRig
from ortho import *
from trees import *
import nodeconfig
from edgeoverlaps import removeEdgeOverlaps
from edgeoverlaps import removeEdgeCrossings
from treeplacement import reattachTrees2
import nearalign
from chainconfig import Chain
from ortho import Compass
from logging import Logger
from logging import LogLevel

class HolaConfig:
    def __init__(self):
        self.DEFAULT_TREE_DIREC = Compass.SOUTH
        """
        Ideal edge length will be a multiple of the average node dimension.
        Set the multiplier here.

        Should be at least two, or else the shape buffer multiplier should not be used.
        """
        self.IEL_MULTIPLIER = 2
        self.ROTATE_FOR_WIDE_ASPECT_RATIO = True
        self.CLASSIC_ACA_FOR_CHAINS = True
        """
        Set the kinds of placements that are favoured for trees.
        """
        self.TREE_PLACEMENT_FAVOUR_CARDINAL = True
        self.TREE_PLACEMENT_FAVOUR_EXTERNAL = True
        self.TREE_PLACEMENT_FAVOUR_ISOLATION = True
        """
        For neighbour stress we scale the ideal edge length.
        """
        self.NBR_STRESS_IEL_SCALAR = 1/20.0
        """
        Options and parameters for the "near alignments" pass:
        """
        self.DO_NEAR_ALIGNMENTS = True
        self.ALIGN_AND_SHAKE_REPS = 2
        self.KINK_WIDTH_SCALAR = 0.5
        self.ALIGNMENT_SCOPE_SCALAR = 2.0
        """
        Options for symmetric tree layout:
        """
        self.RIGID_RANK_SEP = True
        self.TRY_MIRROR_TRIPLES = False
        """
        Logging levels for various stages of the process:
        """
        self.LOG_LEVEL_GENERAL = LogLevel.TIMING
        self.LOG_LEVEL_TREE_PLACEMENT = LogLevel.TIMING

        """
        Do you want node IDs to be labels on the nodes in the final layout?
        """
        self.NODE_IDS_AS_LABELS = False

        """
        Do a final stress reduction with neighbour stress?
        """
        self.DO_FINAL_NEIGHBOUR_STRESS_SHAKE = True

        """
        If three nodes u, v, w in a graph form a triangle -- i.e. a subgraph isomorphic
        to K3 -- then during node configuration we may wish to prevent the "flattening"
        of this triangle, i.e. the configuration of any one of the three nodes in such
        a way that the other two are assigned to opposite compass directions, e.g. assigning
        u and w to be north and south, resp., of v. To prevent this, set this option to
        True.

        One reason to leave this option set to False is e.g. that K4 gets a better (planar)
        layout.
        """
        self.NODE_CONFIG_NO_FLAT_TRIANGLES = False

        """
        To get a dictionary of routing options to be passed to a RoutingRig object,
        call the getRoutingOpts function.

        In many cases it makes sense to think of router parameters as scalar multiples
        of the ideal edge length of the graph. In such cases, you may set the scalars
        in the following dictionary, and you must pass the ideal edge length as iel
        to the getRoutingOpts function.

        If you want to switch some of these off, while
        continuing to use others, simply set their value to None. Where scalars are
        not used we fall back on the defaults dictionary below.
        """
        self.ROUTING_OPT_IEL_SCALARS = {
            'crossingPenalty': 2,
            'shapeBufferDistance': 0.125,
            'segmentPenalty': 0.5
        }
        self.ROUTING_OPT_DEFAULTS = {
            'crossingPenalty': 0,
            'shapeBufferDistance': 8,
            'segmentPenalty': 50
        }
        """
        In, for example, a NORTH-growing tree, an edge between ranks i and i + 1
        will always be allowed to connect only to the south (S) port of a node in
        rank i + 1.

        This setting controls the directions allowed for connection to nodes in
        rank i, as follows:

        0:  only N is allowed
        1:  N, E, W are allowed for the root node if it has exactly one child and
            it is an ordinal placement, otherwise only N
        2:  N, E, W are allowed for all nodes on rank i

        The "CORE" version controls trees attached to a core graph.
        The "PURE" version controls graphs which are themselves trees.
        """
        self.PERMISSIVE_CORE_TREE_ROUTING = 1
        self.PERMISSIVE_PURE_TREE_ROUTING = 2


        """
        How to react if we get a positive water level route with no bends?
        This is indicative of some systematic error, but we may nevertheless
        want to skip over it, and simply mark the path as unusable.
        OR we can even throw all caution to the wind and try to use the path
        anyway.

        The settings are as follows:

        0:  Do not tolerate. Raise an exception and quit immediately.
        1:  Raise an UnusableWaterPath exception, and print a warning.
            This type of exception is caught by a higher level control loop.
        2:  Do not raise any exception, but do print a warning.
        3:  Do not raise any exception, do not print any warning.
        """
        self.ON_POSITIVE_WATER_LEVEL_ROUTE_WITHOUT_BENDS = 1

        """
        Default operation is to evaluate all tree placement options exactly
        by actually carrying out each potential projection sequence, evaluating
        the stress change, and backtracking.
        If speed is favoured over quality, we can instead merely estimate the
        cost of each tree placement. Set this to True if that is desired.
        """
        self.ESTIMATE_TREE_PLACEMENT_COSTS = False
        """
        We can also speed up tree placement by using an estimate of the stress
        costs, in order to choose the primary expansion dimension.
        """
        self.HEURISTIC_CHOICE_FOR_PRIMARY_EXPANSION_DIMENSION = False
        """
        When making heuristic choice of primary expansion dimension, do you want
        to work in the costlier dimension first? (The hope is that the bigger change
        that this represents will already be enough to make room for the tree, and
        you will not have to work in the other dimension at all.)
        """
        self.HCPED_COSTLIER_DIMENSION_FIRST = True

        """
        Ignore all but level zero? Doing so may miss some alternative ways to
        expand a face, but will be faster.
        """
        self.WATER_LEVEL_ZERO_ONLY = False

        """
        Should we use scaling when using stress majorization for neighbour stress layout?
        Recent tests have shown it is faster if we do _not_ use it.
        """
        self.USE_SCALING_IN_MAJORIZATION = False


    def useFastSettings(self, b):
        """
        :param b: boolean
        :return: nothing

        Use all the settings designed for greater speed, if b == True.
        """
        self.ESTIMATE_TREE_PLACEMENT_COSTS = b
        self.HEURISTIC_CHOICE_FOR_PRIMARY_EXPANSION_DIMENSION = b
        self.WATER_LEVEL_ZERO_ONLY = b


    def getKinkWidth(self, avgdim):
        """
        :param avgdim: the average dimension of nodes
        """
        return self.KINK_WIDTH_SCALAR*avgdim

    def getAlignmentScope(self, avgdim):
        return self.ALIGNMENT_SCOPE_SCALAR*avgdim

    def writeConfigCode(self):
        c = "CA" if self.CLASSIC_ACA_FOR_CHAINS else "SB"
        if self.TREE_PLACEMENT_FAVOUR_CARDINAL: c += "_CARD"
        if self.TREE_PLACEMENT_FAVOUR_EXTERNAL: c += "_EXT"
        if self.DO_NEAR_ALIGNMENTS: c += "_NA"
        return c

    def getRoutingOpts(self, iel=0):
        opts = self.ROUTING_OPT_DEFAULTS.copy()
        for p in opts:
            s = self.ROUTING_OPT_IEL_SCALARS.get(p, None)
            if s is not None:
                opts[p] = s * iel
        return opts

    def getShapeBufferDistance(self, iel=0):
        opts = self.getRoutingOpts(iel=iel)
        return opts['shapeBufferDistance']



def hola(G_orig, config=None, logger=None, projLogger=None):
    """
    :param G_orig: the Graph to be laid out
    :param config: optional HolaConfig object
    :param logger: optional Logger
    :param projLogger: optional Logger for tree placement
    :return: nothing

    We perform HOLA layout on the given graph.
    """
    # If no configuration object was passed, we use a default config.
    if config is None:
        config = HolaConfig()
    # If no logger was passed, then we do no logging.
    # Otherwise we set its level to that specified in the config.
    if logger is None:
        logger = Logger('', '')
        logger.level = LogLevel.NONE
    else:
        logger.level = config.LOG_LEVEL_GENERAL
    # You can pass a separate logger for the projections in
    # the tree placement step.
    # Likewise, if you do not pass it then we do no logging for
    # that part of the process; if you do, then we set its level
    # from the config.
    if projLogger is None:
        projLogger = Logger('', '')
        projLogger.level = LogLevel.NONE
    else:
        projLogger.level = config.LOG_LEVEL_TREE_PLACEMENT

    if logger.level >= LogLevel.TIMING:
        logger.startNewTimer('HOLA Layout')

    G_orig.registerMaxID()
    G = G_orig.copy()
    if logger.level >= LogLevel.NODE_IDS_AS_LABELS:
        G.setIDsAsLabels()
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_00_initialGraph", graph=G)

    avgdim = G.computeAvgNodeDim()
    iel = config.IEL_MULTIPLIER * avgdim
    C = prune(G)

    # If it's just a tree, layout and quit.
    if len(C) == 1 and hasattr(C[0], 'rootNode'):
        growthDir = config.DEFAULT_TREE_DIREC
        g = C[0]
        t = Tree(g, g.rootNode)
        t.setHolaConfig(config)
        t.symmetricLayout(growthDir, avgdim/2.0, iel)

        cleanupNodeAppearance(g, config)
        rig = RoutingRig({
            'routing': adg.OrthogonalRouting,
            'nudgeOrthogonalSegmentsConnectedToShapes': True,
            'nudgeSharedPathsWithCommonEndPoint': True
        })
        rig = t.buildOwnRoutingRig(rig, permissive=config.PERMISSIVE_PURE_TREE_ROUTING)
        rig.route(setRoutesInEdges=True)

        t.graph.setPosesInCorrespNodes(G_orig)
        t.graph.setRoutesInCorrespEdges(G_orig)

        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer() # HOLA Layout

        return

    # Otherwise we do have a trunk and trees.
    trunk = C[0]
    treeGraphs = C[1:]

    trunk.setIDsAsLabels()
    trunk.iel = iel

    # We want nodes to be padded through Step 5.
    # The we will remove the padding before routing the edges in Step 6,
    # so that there are "alleyways" between the nodes.
    node_padding = config.getShapeBufferDistance(iel=iel)
    # Save scalar to restore later.
    shapeBufferDistanceScalar = config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance']
    config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance'] = None
    trunk.addPaddingToNodes(node_padding, node_padding)

    # Do an FD layout with NO overlap prevention.
    startT('01_FD', logger)
    trunk.fdlayout(iel, False)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_01_FD", graph=trunk)
    # Do another FD layout WITH overlap prevention.
    startT('02_FDwOP', logger)
    trunk.fdlayout(iel, True)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_02_FDwOP", graph=trunk)

    # Configure nodes in trunk:
    startT('NWACA', logger)
    ccs = adg.CompoundConstraintPtrs()
    aca, rs = nodeconfig.nodewiseACA(trunk, iel, ccs, logger, config)[:2]
    trunk.inferNodeConf()
    trunk.nodeConf.extraGapX = avgdim
    trunk.nodeConf.extraGapY = avgdim
    stopT(logger)

    # Do an FD layout to shake out accumulated stress and regain symmetries.
    startT('04_shakeSE', logger)
    trunk.shakeWithSolidEdges(iel)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_04_shakeSE", graph=trunk)

    # Configure chains in trunk:
    chains = None
    if config.CLASSIC_ACA_FOR_CHAINS:
        startT('Chains_CA', logger)
        aca.createAlignments()
        trunk.moveNodesToRects(rs)
        trunk.inferNodeConf()
        trunk.nodeConf.extraGapX = avgdim
        trunk.nodeConf.extraGapY = avgdim
        trunk.shakeWithSolidEdges(iel)
        stopT(logger)
    else:
        startT('Chains_SB', logger)
        cc = trunk.getChainsAndCycles()
        chains = [Chain(trunk, links) for links in cc['chains']]
        cycles = [Chain(trunk, links, cycle=True) for links in cc['cycles']]
        # We can just put them all together:
        chains.extend(cycles)
        if logger.level >= LogLevel.DEBUG:
            print "Chains in trunk:"
            for chain in chains:
                print chain
                pbs = chain.possibleBendSeqs()
                print '    pbs: %s\n' % pbs
        changes = False
        for chain in chains:
            changes |= chain.takeShapeBasedConfiguration()
        if changes:
            # Project before shaking with solid edges, so that the edges of
            # the chain can be axis aligned first.
            # We do NOT want overlap prevention for the projection, because
            # the new chain configuration constraints may very well reverse
            # one or more orthogonal orderings, and we need them to be free
            # to do that.
            # (Alternatively, could put OP on, but raise the accept level to 1.)
            trunk.project(logger, adg.XDIM, op=False)
            trunk.project(logger, adg.YDIM, op=False)
            trunk.shakeWithSolidEdges(iel)
        stopT(logger)
    # Show result of chain config:
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_05_chainConfig", graph=trunk)

    # Now we remove the node padding.
    trunk.addPaddingToNodes(-node_padding, -node_padding)

    # Route connectors.
    doRouting = True
    if doRouting:
        startT('06_routing', logger)
        # Save shape buffer distance scalar, and then set temporarily to zero.
        sbd_scalar = config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance']
        config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance'] = 0
        # Now do the routing.
        orthoRoute(trunk, logger, config, iel=iel)
        # Restore shape buffer distance.
        config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance'] = sbd_scalar
        stopT(logger)
        if logger.level >= LogLevel.STAGE_GRAPHS:
            logger.writeGML("_06_routing", graph=trunk)

    # Now restore shape buffer distance for routing.
    config.ROUTING_OPT_IEL_SCALARS['shapeBufferDistance'] = shapeBufferDistanceScalar

    # Remove edge overlaps.
    startT('07_EOR', logger)
    Q = removeEdgeOverlaps(trunk)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_07_EOR", graph=Q)

    # Remove edge crossings.
    startT('08_ECR', logger)
    P = removeEdgeCrossings(Q, withConstraints=True)
    P.nodeConf.extraGapX = avgdim
    P.nodeConf.extraGapY = avgdim
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_08_00_ECR", graph=P)

    # Do an FD layout to shake out accumulated stress and regain symmetries.
    startT('08_01_shakeSE', logger)
    P.shakeWithSolidEdges(iel)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_08_01_shakeSE", graph=P)

    # Lay out trees.
    startT('Tree layout', logger)
    trees = []
    growthDir = config.DEFAULT_TREE_DIREC
    for g in treeGraphs:
        g.setIDsAsLabels()
        t = Tree(g, g.rootNode)
        t.setHolaConfig(config)
        trees.append(t)
        t.symmetricLayout(growthDir, avgdim/2.0, iel)
    stopT(logger)

    P.setLogger(projLogger)

    # Compute shortest paths in P so that it can compute stress locally.
    P.computeShortestPaths(iel)

    # Expand for trees.
    if logger.level >= LogLevel.DEBUG:
        print 'Expanding for trees'
    startT('10_tree_placement', logger)
    faceset = reattachTrees2(P, trees, iel, projLogger, config)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_10_00_expanded_for_trees", graph=P)
    # Do an FD layout to shake out accumulated stress and regain symmetries.
    startT('10_01_shakeSE', logger)
    P.shakeWithSolidEdges(iel)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_10_01_shakeSE", graph=P)

    tps = faceset.getAllTreePlacements()

    # Get original graph.
    G = G_orig

    # Give trees the geometry of their chosen placements.
    for tp in tps:
        tp.applyGeometryToTree(iel=iel)

    # Insert the trees into graph P.
    startT('11_0_insertTrees', logger)
    for tp in tps:
        tp.insertTreeIntoGraph(P, iel=iel)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_11_0_insertTrees", graph=P)

    # Shake with neighbour stress.
    nbr_iel = iel * config.NBR_STRESS_IEL_SCALAR
    logger.pushLevel()
    logger.level = LogLevel.DEBUG
    startT('11_5_nbr_stress', logger)
    P.shakeWithSolidEdges(nbr_iel, logger=logger,
                          useNeighbourStress=True, useScaling=config.USE_SCALING_IN_MAJORIZATION)
    stopT(logger)
    logger.popLevel()

    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_11_5_nbr_stress", graph=P)

    # Need list of all buffer nodes, for near-alignments, and for rotation.
    bufferNodes = []
    for tp in tps:
        bufferNodes.extend(tp.bufferNodes)

    if config.DO_NEAR_ALIGNMENTS:
        # Do the final alignments step on P.
        startT('12_0_nearAlignments', logger)
        ig = bufferNodes

        startT('init_alignment_table', logger)
        atab = nearalign.AlignmentTable(P, ignore=ig)
        for tp in tps:
            a = tp.computeAlignedSets()
            for s in a['h']:
                atab.addAlignedNodes(s, nearalign.AlignmentFlags.HALIGN)
            for s in a['v']:
                atab.addAlignedNodes(s, nearalign.AlignmentFlags.VALIGN)
        stopT(logger)

        for rep in range(config.ALIGN_AND_SHAKE_REPS):
            startT('do_near_alignments', logger)
            nearalign.doNearAlignments(avgdim, logger, P, atab, config, ignore=ig)
            stopT(logger)
            if logger.level >= LogLevel.STAGE_GRAPHS:
                logger.writeGML("_12_0_nearAlignments_%02d" % rep, graph=P)

            # Shake with neighbour stress.
            startT('12_5_nbr_stress', logger)
            P.shakeWithSolidEdges(nbr_iel, logger=logger,
                                  useNeighbourStress=True, useScaling=config.USE_SCALING_IN_MAJORIZATION)
            stopT(logger)
            if logger.level >= LogLevel.STAGE_GRAPHS:
                logger.writeGML("_12_5_nbr_stress_%02d" % rep, graph=P)
        stopT(logger)

    # Rotate and translate
    startT('13_rotate_and_translate', logger)
    rotateAndTranslate(P, faceset, tps, bufferNodes, iel, nbr_iel, config, logger=logger, nbr=True)
    stopT(logger)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_13_9_rotateAndTranslate", graph=P)

    # Set node positions.
    P.setPosesInCorrespNodes(G)

    # Clean up node appearance.
    startT('Final cleanup and routing', logger)
    cleanupNodeAppearance(G, config)

    # Add aesthetic bends in chains.
    # Subtle point: chains' BPs belong to trunk graph, so must
    # update their positions from P.
    P.setPosesInCorrespNodes(trunk)

    # Now can set the route points.
    G.clearAllRoutes()
    if chains is not None:
        for chain in chains:
            chain.addRoutePointsInGraph(G)
        G.buildRoutes()

    # Set up a router for all other edges.
    rig = trunk.buildRigForEdgesNeedingRoute()
    opts = config.getRoutingOpts(iel=iel)
    rig.setOpts(opts)
    for tree in trees:
        tree.addInternalNodesAndEdgesToRoutingRig(rig, core=G,
                                    permissive=config.PERMISSIVE_CORE_TREE_ROUTING)
    startT('Routing', logger)
    rig.route(setRoutesInEdges=True)
    stopT(logger)
    rig.router.outputInstanceToSVG('testOut/cleanup/routes')
    for tree in trees:
        tree.setEdgeRoutesInGraph(G)

    stopT(logger) # Final cleanup and routing

    G.setPosesInCorrespNodes(G_orig)
    G.setRoutesInCorrespEdges(G_orig)

    if logger.level >= LogLevel.TIMING:
        logger.stopLastTimer() # HOLALayout

    # END HOLA ################################################################


def startT(name, logger):
    if logger.level >= LogLevel.TIMING:
        logger.startNewTimer(name)

def stopT(logger):
    if logger.level >= LogLevel.TIMING:
        logger.stopLastTimer()


def rotateAndTranslate(G, faceset, treePlacements, bufferNodes, iel_for_constraints, iel_for_shake, config, logger=None, nbr=False):
    x, y, w, h = G.boundingBoxXYWH(ignore=bufferNodes)
    # Rotate 90 degrees if necessary so the layout is not taller than it is wide.
    if w < h and config.ROTATE_FOR_WIDE_ASPECT_RATIO:
        nw = faceset.numTreesWithGrowthDirec(Compass.WEST, scaleBySize=True)
        ne = faceset.numTreesWithGrowthDirec(Compass.EAST, scaleBySize=True)
        # Choose direction so that there are not more north-growing trees than south.
        if nw > ne:
            #print 'ROTATE 90 ACW'
            G.rotate90acw(iel_for_constraints, iel_for_shake, logger=logger, nbr=nbr)
        else:
            #print 'ROTATE 90 CW'
            G.rotate90cw(iel_for_constraints, iel_for_shake, logger=logger, nbr=nbr)
    # Otherwise rotate 180 degrees if necessary so there are not more
    # north-growing trees than south.
    else:
        nn = faceset.numTreesWithGrowthDirec(Compass.NORTH, scaleBySize=True)
        ns = faceset.numTreesWithGrowthDirec(Compass.SOUTH, scaleBySize=True)
        if nn > ns:
            #print 'ROTATE 180'
            G.rotate180(iel_for_constraints, iel_for_shake)
    # Put upper-left corner of layout bounding box at (0, 0).
    x, y, w, h = G.boundingBoxXYWH()
    G.translateBy((-x, -y))


def cleanupNodeAppearance(G, config):
    G.doUseYEdColour()
    if config.NODE_IDS_AS_LABELS:
        G.setIDsAsLabels()
    else:
        G.eraseAllLabels()

def orthoRoute(G, logger, config, iel=0):
    """
    Do an orthogonal routing of the edges in a graph having no leaves
    :param G: a graph with no leaves
    :return: nothing
    """
    router = adg.Router(adg.OrthogonalRouting)
    router.setRoutingOption(adg.nudgeSharedPathsWithCommonEndPoint, False)
    rig = RoutingRig({})
    opts = config.getRoutingOpts(iel=iel)
    rig.setOpts(opts, router=router)
    connDirsAllowed = {}
    edgeReptoConnRef = {}
    maxNodeID = -1
    for ID in G.nodes:
        node = G.nodes[ID]
        poly = node.libavoidPoly()
        sr = adg.ShapeRef(router, poly, ID)
        maxNodeID = max(maxNodeID, ID)
    baseEdgeID = maxNodeID + 1
    for i, rep in enumerate(G.edges.keys()):
        edge = G.edges[rep]
        ID = baseEdgeID + i
        cr = adg.ConnRef(router, ID)
        edgeReptoConnRef[rep] = cr
        src, tgt = edge.libavoidConnEnds()
        cr.setEndpoints(src, tgt)
        connDirsAllowed[rep] = {
            edge.srcID: adg.ConnDirAll,
            edge.tgtID: adg.ConnDirAll
        }
    # We will need a function to tell us the direction in which a routed
    # connector departs a node.
    def departureDirec(e, u):
        """
        e: an Edge
        u: a Node being one endpt of the edge
        We return the Compass direction in which the routed connector for edge e
        departs node u.
        """
        cr = edgeReptoConnRef[repr(e)]
        polyline = cr.displayRoute()
        if u.ID == e.srcID:
            pts = polyline.ps
        else:
            assert(u.ID == e.tgtID)
            pts = list(reversed(polyline.ps))
        firstpair = [
            (p.x, p.y) for p in pts[:2]
        ]
        p0, p1 = firstpair
        d = Compass.cardinalDirection(p0, p1)
        return d

    # Now route. We may need to route multiple times to ensure that at least /two/
    # sides of each node are being used, but in theory we should never have to route
    # more than /five/ times.
    # To test that theory we use an infinite loop with counter and assertion, instead
    # of a mere for-loop which would fail silently.
    numRoutings = 0
    while True:
        #print 'numRoutings: %s' % numRoutings
        router.processTransaction()
        router.outputInstanceToSVG("output/routing%d" % numRoutings)
        numRoutings += 1
        assert(numRoutings <= 5)
        # Are there any nodes having all of their edges routed
        # out of just one side?
        # We need at least two sides of each node to be used, or else
        # our planarisation will have leaves.
        pseudoLeaves = []
        for node in G.nodes.values():
            edges = node.edges.values()
            # Since we assume there are no actual leaves in the graph:
            assert(len(edges) >= 2)
            e0 = edges[0]
            d0 = departureDirec(e0, node)
            for e in edges[1:]:
                d = departureDirec(e, node)
                if d != d0: break
            else:
                pseudoLeaves.append(node)
        if len(pseudoLeaves) == 0:
            break
        else:
            for node in pseudoLeaves:
                edges = node.edges.values()
                d0 = departureDirec(edges[0], node)
                # Choose an edge whose connection directions we will restrict.
                # We must choose only from among those edges that are still permitted
                # at least two directions, since we will forbid one of these directions,
                # and at least one must remain!
                loneDirecs = Compass.libavoidVisibility.values()
                edges = filter(
                    lambda e: connDirsAllowed[repr(e)][node.ID] not in loneDirecs,
                    edges
                )
                # Choose edge to first neighbour which is /not/ in direc d0 from the node,
                # if possible; else just accept the last one.
                for e0 in edges:
                    v = e0.otherEnd(node)
                    d = Compass.cardinalDirection(node, v)
                    if d != d0: break
                # Get the connector.
                rep = repr(e0)
                cr = edgeReptoConnRef[rep]
                # Start with directions allowed last time:
                available = connDirsAllowed[rep][node.ID]
                # XOR with d0, so it is no longer available:
                available ^= Compass.libavoidVisibility[d0]
                connDirsAllowed[rep][node.ID] = available
                # Set the new connend:
                pt = adg.Point(node.x, node.y)
                ce = adg.ConnEnd(pt, available)
                if node.ID == e0.srcID:
                    cr.setSourceEndpoint(ce)
                else:
                    assert(node.ID == e0.tgtID)
                    cr.setDestEndpoint(ce)

    # Finally, store the connector routes in the Edge objects.
    for rep in edgeReptoConnRef:
        cr = edgeReptoConnRef[rep]
        polyline = cr.displayRoute()
        if logger.level >= LogLevel.DEBUG:
            print 'Edge: %s' % rep
            print ' '*4 + ', '.join([
                '(%.2f, %.2f)' % (p.x, p.y) for p in polyline.ps
            ])
        # On some examples we were getting a three-segment
        # route in which the middle vertical segment was about 10e-15 in length.
        # We remove such segments here.
        EPSILON = 0.001
        pts = [(p.x, p.y) for p in polyline.ps]
        interior = pts[1:-1]
        N, i = len(interior), 0
        pts1 = []
        def samePts(p, q):
            px, py = p
            qx, qy = q
            return abs(px - qx) < EPSILON and abs(py - qy) < EPSILON
        while i + 1 < N:
            p, q = interior[i:i+2]
            if samePts(p, q):
                i += 2
            else:
                pts1.append(p)
                i += 1
        if i < N:
            pts1.append(interior[i])
        pts = pts[:1] + pts1 + pts[-1:]
        # Store the route in the Edge.
        edge = G.edges[rep]
        edge.route = pts

