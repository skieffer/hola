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
import assignments.assign as assign

from logging import LogLevel
from collections import deque as Deque


def nodewiseACA(G, L, ccs, logger, config):
    """
    :param G: graph to which nodewise ACA is to be applied
    :param L: ideal edge length
    :param ccs: compound constraints vector containing any existing constraints,
                and to which new ones should be added
    :return: the ACA object that we build here
    """
    # Get CoLa structures for G.
    rs, es, ix = G.writeRsEsIx()

    if logger.level >= LogLevel.DEBUG:
        nodes = G.nodes.values()
        nodes.sort(key=lambda u: u.ID)
        print 'ix map for NWACA:'
        for u in nodes:
            print '    ID %s --> index %s' % (u.ID, ix(u))

    # Set node labels for debugging purposes.
    if logger.level >= LogLevel.NODE_IDS_AS_LABELS:
        G.setIDsAsLabels()

    # We'll need to do an FD layout once in a while.
    def fdlayout(op=True):
        alg = adg.ConstrainedFDLayout(rs, es, L, op)
        alg.setConstraints(ccs)
        alg.run()
        G.moveNodesToRects(rs)

    # Sort nodes by degree.
    nodes = G.nodes.values()
    nodes.sort(key=lambda u: u.degree, reverse=True)
    # Exclude all degree-2 nodes.
    for i, u in enumerate(nodes):
        if u.degree == 2: break
    else:
        i += 1
    nodes = nodes[:i]
    # Create an ACA layout object.
    aca = adg.ACALayout3(rs, es, ccs, L, True)
    # Go through nodes one by one, in order, giving each one the best possible
    # configuration that is feasible with the existing constraints.
    if logger.level >= LogLevel.FINER_STAGE_GRAPHS:
        numSteps = 0
    ptr = 0
    mightNeedToShake = True  # says whether we can reduce stress and try a node again
    while ptr < len(nodes):
        centre = nodes[ptr]
        if logger.level >= LogLevel.DEBUG:
            print 'Node: %d, index = %d, deg = %d' % (
                centre.ID, ix(centre), centre.degree
            )
        # Compute the possible assignments of neighbours to compass directions,
        # in order of descending desirability.
        # First make sure the Node objects have the most up-to-date positions.
        G.moveNodesToRects(rs)
        # Now compute list of possible assignments.
        asgns = Deque(assign.getAssignmentsForNode(centre))
        # Filter out assignments that would result in flat triangles, if we
        # want to disallow those.
        if config.NODE_CONFIG_NO_FLAT_TRIANGLES:
            # Precompute the adjacency matrix.
            am = G.getAdjacencyMatrix()
            # Write filter function.
            def makesFlatTriangle(asgn):
                for i in range(2):
                    s0, s1 = asgn.semis[i], asgn.semis[i + 2]
                    for nbr0 in s0:
                        ID0 = nbr0.ID
                        for nbr1 in s1:
                            ID1 = nbr1.ID
                            if am[ID0].get(ID1, False):
                                return True
                return False
            asgns = Deque(filter(lambda a: not makesFlatTriangle(a), asgns))
        # Logging
        if logger.level >= LogLevel.DEBUG:
            print '%d possible assignments' % len(asgns)
        # If no assignments, move on.
        if len(asgns) == 0:
            ptr += 1
            continue
        # Start trying them, and quit either when one works or when we run out.
        success = False
        oas = adg.OrderedAlignmentPtrs()
        while asgns:
            asgn = asgns.popleft()
            if logger.level >= LogLevel.DEBUG:
                print '    %s' % asgn
            oas.clear()
            for i, semi in enumerate(asgn.semis):
                # If the assignment places no node on semiaxis i, continue
                # to next semiaxis.
                if len(semi) == 0: continue
                for nbr in semi:
                    ci, ni = (ix(centre.ID), ix(nbr.ID))
                    # If these nodes are already logically aligned, then skip this OA.
                    # (This is not just to save time. In fact if you give VPSC redundant
                    #  equality constraints it will mark the second one as unsatisfiable.
                    #  This is because once one constraint is active, it will think the
                    #  other is still inactive and violated. To satisfy that one it
                    #  will try to split the block to which the two variables already
                    #  belong, and fail because there are only equality constraints
                    #  along the path between them. Since we are calling ACA's
                    #  'allOrNothing' method, this one "unsatisfiable" constraint will
                    #  cause the entire node arrangement to fail.)
                    if aca.nodesAreAligned(ci, ni): continue
                    # Else grab the appropriate ACASEPFLAG and create the OA.
                    flag = [adg.ACAEAST, adg.ACASOUTH, adg.ACAWEST, adg.ACANORTH][i]
                    oa = aca.initOrdAlign(ci, ni, flag, getEdgeIndex(ci, ni, es))
                    oas.push_back(oa)
            if logger.level >= LogLevel.DEBUG:
                data = (centre.ID,) + asgn.prettyIDTuple()
                name = 'before_%s.%s.%s.%s.%s' % data
                writeAONcppTest(name, rs, es, ccs, L, oas)
                G.moveNodesToRects(rs)
                gml = G.writeGML()
                f = open('testOut/nwaca_debug/%s.gml'%name,'w')
                f.write(gml)
                f.close()
            success = aca.applyOAsAllOrNothing(oas)
            if logger.level >= LogLevel.DEBUG:
                print '    fail' if not success else '    success'
            if success: break
        # Now we have either found an assignment that works, or tried them all and
        # none of them worked.
        if success:
            # If any assignment was successful, the ACA object automatically adds the
            # new compound constraints to the ccs vector, which was passed by reference,
            # so we do not need to update that ourselves.
            ptr += 1
            mightNeedToShake = True
            if logger.level >= LogLevel.FINER_STAGE_GRAPHS:
                G.moveNodesToRects(rs)
                logger.writeGML("_03_NWACA_%03d_config_node_%d" % (numSteps, centre.ID), graph=G)
                numSteps += 1
        elif not success and mightNeedToShake:
            # If we were not able to configure this node, it may be that relieving
            # stress in the graph will permit us to configure it.
            # So run an FD layout and try this node again.
            fdlayout(op=True)
            mightNeedToShake = False
            if logger.level >= LogLevel.FINER_STAGE_GRAPHS:
                logger.writeGML("_03_NWACA_%03d_shake"%numSteps, graph=G)
                numSteps += 1
        else:
            # If we have already tried relieving stress once, and we /still/ couldn't
            # apply any assignment, then we give up on this node and move on to the next.
            # Note that we must leave mightNeedToShake equal to False, because since we have
            # tried shaking once already, there is no reason to shake again until at least
            # one more node has been configured.
            ptr += 1
    # Accept final positions.
    G.moveNodesToRects(rs)
    if logger.level >= LogLevel.STAGE_GRAPHS:
        logger.writeGML("_03_NWACA_done", graph=G)
    return (aca, rs, es)

def getEdgeIndex(r1, r2, es):
    """
    Find the index of an edge having given endpoints.
    :param r1: index into rs array
    :param r2: index into rs array
    :param es: array of cola edges
    :return: the index of the (first) edge in es having r1 and r2 as its
             endpoints, in either order; or -1 if no edge is found
    """
    # Ensure that r1 <= r2.
    if r2 < r1: r1, r2 = (r2, r1)
    for i, e in enumerate(es):
        if r1 == min(e) and r2 == max(e):
            break
    else:
        i = -1
    return i
