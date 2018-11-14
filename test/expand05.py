
"""
14 Nov 2018
This test exposes a flaw in the face expansion code, which has been fixed in the
new C++ version in Adaptagrams.
"""

import sys, os
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

from hola.gmlparse import buildGraph
import hola.logging as logging
from hola.hola import HolaConfig
from hola.treePruner3 import prune
from hola.trees import *
from hola.edgeoverlaps import removeEdgeOverlaps
from hola.edgeoverlaps import removeEdgeCrossings
from hola.faces import FaceSet
from hola.faces import TreePlacement
from hola.ortho import Compass

def main():

    # Build the graph.
    name = 'special/expand05'
    with open('graphs/%s.gml' % name) as f:
        gml = f.read()
    G = buildGraph(gml)
    G.registerMaxID()

    # Compute IEL.
    avgdim = G.computeAvgNodeDim()
    IEL = 2 * avgdim

    # Set up logging.
    groupFolder = 'output'
    parts = name.split('/')
    uniqueName = parts[-1]
    if len(parts) == 2: groupFolder += '/' + parts[0]
    final = groupFolder + '/final'
    if not os.path.exists(final):
        os.makedirs(final)
    logger = logging.Logger('output/'+name, uniqueName)
    projLogger = logging.Logger('output/'+name+'/proj', uniqueName)

    # Get a HolaConfig.
    config = HolaConfig()

    # Peel
    C = prune(G)
    core = C[0]
    core.setIDsAsLabels()
    core.iel = IEL

    treeGraphs = C[1:]
    t0 = None
    growthDir = config.DEFAULT_TREE_DIREC
    for g in treeGraphs:
        g.setIDsAsLabels()
        t = Tree(g, g.rootNode)
        t.setHolaConfig(config)
        t.symmetricLayout(growthDir, avgdim/2.0, IEL)
        rid = g.rootNode.ID
        if rid == 0: t0 = t

    # Planarise the core.
    Q = removeEdgeOverlaps(core)
    P = removeEdgeCrossings(Q, withConstraints=True)
    logger.writeGML("_00_planar_layout", graph=P)

    # Compute the faces of the core.
    faceset = FaceSet(P, logger, config)
    # List all possible placements for the tree.
    tps0 = faceset.listAllPossibleTreePlacements(t0)
    # Find the one we want.
    tp0 = None
    for tp in tps0:
        if tp.placementDirec == Compass.SE and tp.growthDirec == Compass.EAST:
            tp0 = tp
            break

    assert isinstance(tp0, TreePlacement)
    # Now try to expand for the placement at Node 0.
    ps = tp0.getBestProjSeq(IEL)
    P.applyProjSeq(ps, logger, opx=True, opy=True, solidEdgesX=True, solidEdgesY=True)


if __name__ == '__main__':
    main()
    #import cProfile
    #cProfile.run('main()')

