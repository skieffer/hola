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

from faces import FaceSet
from ortho import Compass
from graphs import Node
from graphs import Edge
from logging import LogLevel

def reattachTrees2(trunk, trees, iel, logger, config):
    """
    :param trunk: a Graph object, being the trunk of a graph, with an existing layout
    :param trees: a list of Tree objects, being the trees of the graph, with existing layout
    :param iel: ideal edge length for the graph
    :return: the FaceSet for the trunk, in which the Faces will have all the TreePlacements
             representing the placements that we make

    We add to trunk a node representing each tree, having the tree's bounding
    box dimensions. Each such node contains a reference to the tree that it
    represents, in the field, 'actualTree'.
    """
    # Compute the faces of the trunk.
    faceset = FaceSet(trunk, logger, config)
    # Sort the trees into the order in which we want to reattach them.
    # For now we try placing those with biggest perimeter first.
    trees.sort(key=lambda t: t.graph.bboxPerimeter(), reverse=True)
    # Reattach the trees one by one.
    if logger.level >= LogLevel.TIMING:
        logger.startNewTimer('trees')
    if logger.level >= LogLevel.PROGRESS:
        placecon = 'Placements considered:\n'
    for tree in trees:
        if logger.level >= LogLevel.TIMING:
            logger.startNewTimer('tree')
        # List all possible placements.
        tps = faceset.listAllPossibleTreePlacements(tree)
        # Choose the best one.
        tp = chooseBestPlacement(tps, iel, logger, config,
                                 favourCardinal=config.TREE_PLACEMENT_FAVOUR_CARDINAL,
                                 favourExternal=config.TREE_PLACEMENT_FAVOUR_EXTERNAL,
                                 favourIsolation=config.TREE_PLACEMENT_FAVOUR_ISOLATION
        )
        # Project, making room for the treenode.
        ps = tp.getBestProjSeq(iel)
        if logger.level >= LogLevel.PROGRESS:
            placecon += '\n'
            for tp1 in tps:
                placecon += '    %s\n' % repr(tp1)
        if logger.level >= LogLevel.TIMING:
            logger.startNewTimer('projection')
        trunk.applyProjSeq(ps, logger, opx=True, opy=True, solidEdgesX=True, solidEdgesY=True)
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
        # Insert treenode.
        tp.insertTreeNode(iel)
        if logger.level >= LogLevel.TIMING:
            logger.stopLastTimer()
    if logger.level >= LogLevel.TIMING:
        logger.stopLastTimer()
    if logger.level >= LogLevel.PROGRESS:
        print placecon
    if logger.level >= LogLevel.TIMING and False:
        td = logger.rootTimer.getTimingDataRec()
        np = logger.numProjections
        print '%d projections, total' % np
        print 'TOTAL TIMES:'
        print td.write()
        nt = len(trees)
        print '%d trees' % nt
        if nt > 0:
            print 'avg %.2f projections per tree' % (float(np)/float(nt))
            print 'AVG TIMES:'
            print td.write(div=nt)
    return faceset


def chooseBestPlacement(tps, iel, logger, config,
                        favourCardinal=True, favourExternal=True, favourIsolation=True):
    """
    :param tps: list of TreePlacement objects
    :param iel: ideal edge length for the graph
    :return: the best one
    """
    bestPlacement = None

    def cmpCardinal(p, q):
        pc = p.placementDirec in Compass.cwCards
        qc = q.placementDirec in Compass.cwCards
        if pc and not qc:
            return -1
        elif qc and not pc:
            return 1
        else:
            return 0

    def cmpExternal(p, q):
        pe = p.face.external
        qe = q.face.external
        if pe and not qe:
            return -1
        elif qe and not pe:
            return 1
        else:
            return 0

    def cmpIsolation(p, q):
        return p.getNumPotentialNbrs() - q.getNumPotentialNbrs()

    if favourCardinal:
        tps.sort(cmp=cmpCardinal)
        # How many of the placements are in a cardinal direction?
        # Due to sorting, these all come first, if any.
        numCardinal = 0
        for tp in tps:
            if tp.placementDirec in Compass.cwCards:
                numCardinal += 1
            else:
                break
        if numCardinal == 1:
            # There is a unique cardinal placement. Choose it.
            bestPlacement = tps[0]
        else:
            # If there are several cardinal placements, then we choose only from among them.
            if numCardinal > 1:
                tps = tps[:numCardinal]

    if bestPlacement is None and favourExternal:
        tps.sort(cmp=cmpExternal)
        # Consider how many placements are in the external face.
        numExternal = 0
        for tp in tps:
            if tp.face.external:
                numExternal += 1
            else:
                break
        if numExternal == 1:
            # There is a unique external placement. Choose it.
            bestPlacement = tps[0]
        else:
            # If there are several external placements, then we choose only from among them.
            if numExternal > 1:
                tps = tps[:numExternal]

    if bestPlacement is None and favourIsolation:
        # Sort tps by number of potential nbrs.
        nbrnums = [(tp, tp.getNumPotentialNbrs()) for tp in tps]
        nbrnums.sort(key=lambda p: p[1])
        # Get all those that have minimal number.
        m = nbrnums[0][1]
        i = 0
        while i < len(tps) and nbrnums[i][1] == m: i += 1
        minimal = [p[0] for p in nbrnums[:i]]
        numMinimal = i
        if numMinimal == 1:
            # There is a unique placement with minimal number of potential neighbours. Choose it.
            bestPlacement = minimal[0]
        else:
            # We choose only from among the placements with minimal number of potential nbrs.
            tps = minimal

    if bestPlacement is None:
        # Finally, we come to the case in which we must evaluate the cost of each remaining
        # potential placement, and choose the cheapest one.
        for tp in tps:
            if config.ESTIMATE_TREE_PLACEMENT_COSTS:
                tp.estimateCost(iel, logger, use_old_heuristic=config.USE_OLD_COST_ESTIMATE_HEURISTIC)
            else:
                tp.evaluateCost(iel)
        bestPlacement = min(tps, key=lambda tp: tp.cost)

    assert(bestPlacement is not None)
    return bestPlacement
