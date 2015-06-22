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
from ortho import Compass

class NodeConfig:
    """
    Maybe this class can make it easier to manage the constraints
    on Nodes than the current PyConstraints class...?
    """

    def __init__(self, G):
        """
        :param G: the Graph for which we represent an orthogonal configuration
        """
        self.graph = G
        # We initialise a sparse data structure to record configured compass
        # directions from one node to another.
        # It is sparse not because we want to save space but because we don't want
        # to have a fixed, rigid number of nodes.
        # Should not be queried directly, but instead via methods defined below.
        # Format: smaller node ID points to dict in which larger node ID points to
        #         compass direction.
        self.direcs = {}
        # An extra gap to be added to all sepco gaps:
        self.extraGapX = 0
        self.extraGapY = 0

    def rotate90cw(self):
        for a in self.direcs.values():
            for ID, t in a.items():
                d = t[0]
                q, r = d/8, d%8
                rp = Compass.cw90(r)
                dp = 8*q + rp
                a[ID] = (dp, t[1])

    def rotate90acw(self):
        for a in self.direcs.values():
            for ID, t in a.items():
                d = t[0]
                q, r = d/8, d%8
                rp = Compass.acw90(r)
                dp = 8*q + rp
                a[ID] = (dp, t[1])

    def rotate180(self):
        for a in self.direcs.values():
            for ID, t in a.items():
                d = t[0]
                q, r = d/8, d%8
                rp = Compass.flip(r)
                dp = 8*q + rp
                a[ID] = (dp, t[1])

    def removeNode(self, u):
        """
        :param u: a Node
        :return: nothing

        We remove all records involving the node.
        """
        d = self.direcs
        if d.has_key(u.ID):
            del d[u.ID]
        for e in d.values():
            if e.has_key(u.ID):
                del e[u.ID]

    def copy(self, H):
        """
        :param H: the Graph with which you want the copy to be associated
        :return: a copy of this object
        """
        nc = NodeConfig(H)
        nc.extraGapX = self.extraGapX
        nc.extraGapY = self.extraGapY
        for id1 in self.direcs:
            a = self.direcs[id1]
            b = {}
            for id2 in a:
                b[id2] = a[id2]
            nc.direcs[id1] = b
        return nc

    def setDirec(self, u, v, d, alignOffsets=(0,0)):
        """
        :param u: a Node object
        :param v: a Node object
        :param d: a Compass direction. You may add 8 to signal that the sepco should be exact.
        :return: nothing

        We record that the direction from u to v should be d.
        """
        q, r = d/8, d%8
        id1, id2 = min(u.ID, v.ID), max(u.ID, v.ID)
        if id1 == v.ID:
            r = Compass.flip(r)
            alignOffsets = tuple(reversed(alignOffsets))
        a = self.direcs.get(id1, {})
        a[id2] = (r+8*q, alignOffsets)
        self.direcs[id1] = a

    def getDirec(self, u, v):
        """
        :param u: a Node object
        :param v: a Node object
        :return: the configured Compass direction from u to v if any, else None
        """
        d = None
        id1, id2 = min(u.ID, v.ID), max(u.ID, v.ID)
        a = self.direcs.get(id1, None)
        t = None
        if a is not None:
            t = a.get(id2, None)
        if t is not None:
            d, alignOffsets = t
            q, r = d/8, d%8
            if id1 == v.ID: r = Compass.flip(r)
            d = r+8*q
        return d

    def free(self, u, v):
        """
        Remove any existing configuration of nodes u and v.
        :param u: a Node object
        :param v: a Node object
        :return: the Compass direction that was removed, or None if there was none
        """
        d = self.getDirec(u, v)
        if d is not None:
            id1, id2 = min(u.ID, v.ID), max(u.ID, v.ID)
            a = self.direcs[id1]
            del a[id2]
        return d

    def buildPCS(self):
        """
        Build and return Python constraint objects
        :return: ordered pair (acs, scs) being a list of AlignCo objects and
                 a list of SepCo objects
        """
        acs = []
        scs = []
        for id1 in self.direcs.keys():
            a = self.direcs[id1]
            for id2 in a.keys():
                # Get the desired compass direction d.
                t = a[id2]
                d, alignOffsets = t
                q, r = d/8, d%8
                # Get nodes.
                u1, u2 = self.graph.nodes[id1], self.graph.nodes[id2]
                # Make list of directions.
                direcs = [r] if r in Compass.cwCards else Compass.components[r]
                for d in direcs:
                    # Determine which node is left and which is right.
                    if d in Compass.increasing:
                        ul, ur = u1, u2
                    else:
                        ul, ur = u2, u1
                    # Determine the dimensions for the alignment and separation,
                    # and the gap for the separation.
                    if d in Compass.horizontal:
                        adim, sdim = adg.YDIM, adg.XDIM
                        extra = self.extraGapX if q == 0 else 0
                        gap = (u1.w + u2.w) / 2.0 + extra
                    else:
                        adim, sdim = adg.XDIM, adg.YDIM
                        extra = self.extraGapY if q == 0 else 0
                        gap = (u1.h + u2.h) / 2.0 + extra
                    # Decide exactness.
                    exact = (q == 1)
                    # Create alignment only if r is a cardinal direction.
                    if r in Compass.cwCards:
                        # Create and append the alignment constraint.
                        acs.append(AlignCo(adim, shapes=[
                            [u1, alignOffsets[0]], [u2, alignOffsets[1]]
                        ]))
                    # Create separation constraint.
                    scs.append(SepCo(sdim, ul, ur, gap, exact=exact))
        return (acs, scs)

    def buildCCS(self, ix=None):
        """
        :param ix: optional callable to be applied to node IDs
        :return: a Cola constraint objects vector
        """
        acs, scs = self.buildPCS()
        ccs = adg.CompoundConstraintPtrs()
        for pc in acs + scs:
            pcCCs = pc.buildCCs(ix=ix)
            for cc in pcCCs:
                ccs.push_back(cc)
        return ccs

    def writeCpp(self, ix=None):
        s = ''
        acs, scs = self.buildPCS()
        if acs or scs:
            s += 'CompoundConstraints ccs;\n'
            if acs:
                s += 'AlignmentConstraint *algn = NULL;\n'
            if scs:
                s += 'SeparationConstraint *sep = NULL;\n'
            for c in acs:
                s += c.writeCpp(ix=ix)
                s += 'ccs.push_back(algn);\n'
            for c in scs:
                s += c.writeCpp(ix=ix)
                s += 'ccs.push_back(sep);\n'
        return s

    # FIXME: Instead of passing nextIDs back and forth, should be using an ID Dispenser.
    def writeDunnartSVG(self, nextID):
        """
        :param nextID: everything in Dunnart SVG has an ID; this says what the
                       next available one is
        :return: SVG lines for the constraints implied by this node configuration object
        """
        acs, scs = self.buildPCS()
        s = ''
        for ac in acs:
            t, nextID = ac.writeDunnartSVG(nextID)
            s += t
        for sc in scs:
            t, nextID = sc.writeDunnartSVG(nextID)
            s += t
        return s

class OrdAlign:

    def __init__(self, left, right, sepdir):
        self.left = left
        self.right = right
        self.sepdir = sepdir

    def buildOA(self, aca):
        oa = aca.initOrdAlign(self.left, self.right, self.sepdir)
        return oa

    def writeCpp(self, name='oa'):
        s = ''
        s += '%s = aca->initOrdAlign(%d, %d, %s);\n' % (
            name, self.left, self.right, {
                adg.ACASOUTH: 'ACASOUTH',
                adg.ACAEAST: 'ACAEAST'
            }[self.sepdir]
        )
        return s

class AlignCo:

    def __init__(self, dim, shapes=[], intcoords=False):
        self.dim = dim
        if intcoords:
            self.shapeOffsetPairs = [
                (u, int(round(o))) for u, o in shapes
            ]
        else:
            self.shapeOffsetPairs = shapes

    def __repr__(self):
        s = ''
        d = {adg.XDIM: 'x', adg.YDIM: 'y'}[self.dim]
        s += 'Align %s:' % {adg.XDIM: '|', adg.YDIM: '--'}[self.dim]
        for u, o in self.shapeOffsetPairs:
            s += ' (%d, %.2f)' % (u.ID, o)
        return s

    def rotateCW(self, n):
        """
        :param n: an integer in the set {1, 2, 3}
        :return: nothing

        We rotate the constraint 90*n degrees clockwise.
        """
        assert(1 <= n <= 3)
        dim = self.dim
        if (dim == adg.XDIM and n == 3) or (dim == adg.YDIM and n == 1):
            # These are the two cases in which all offsets must be negated.
            self.shapeOffsetPairs = [
                (u, -o) for u, o in self.shapeOffsetPairs
            ]
        if n != 2:
            self.dim = 1 - dim

    def buildCCs(self, rs=None, ix=None):
        algn = adg.AlignmentConstraint(self.dim)
        for u, off in self.shapeOffsetPairs:
            i = u.ID if ix is None else ix(u.ID)
            algn.addShape(i, off)
        return [algn]

    def writeCpp(self, name='algn', ix=None):
        s = ''
        s += '%s = new AlignmentConstraint(%s);\n' % (name, {
            adg.XDIM: 'XDIM',
            adg.YDIM: 'YDIM'
        }[self.dim])
        for u, off in self.shapeOffsetPairs:
            i = u.ID if ix is None else ix(u)
            s += '%s->addShape(%d, %.2f);\n'%(
                name, i, off
            )
        return s

    def writeDunnartSVG(self, nextID):
        s = ''
        if self.dim == adg.YDIM:
            # "horizontal alignment"
            direc = 101
            coord = lambda u: u.y
        else:
            # "vertical alignment"
            direc = 100
            coord = lambda u: u.x
        pos = sum([
            coord(u) - off
            for u, off in self.shapeOffsetPairs
        ]) / float(len(self.shapeOffsetPairs))
        GID = nextID
        s += '<dunnart:node type="guideline" position="%.2f" direction="%d" id="%d">\n'%(
            pos, direc, GID
        )
        nextID += 1
        ap = {adg.YDIM: 1, adg.XDIM: 4}[self.dim]
        for u, off in self.shapeOffsetPairs:
            # TODO: take acct of offset?
            s += '    <dunnart:node type="constraint" relType="alignment" '
            s += 'alignmentPos="%d" objOneID="%d" constraintID="%d"/>\n' % (
                ap, u.ID, GID
            )
        s += '</dunnart:node>\n'
        return (s, nextID)

class SepCo:

    def __init__(self, dim, left, right, gap, exact=False, intcoords=False):
        """
        :param dim: adg.XDIM or adg.YDIM
        :param left: a Node
        :param right: a Node
        :param gap: float
        :param exact: bool
        """
        self.dim = dim
        self.left = left
        self.right = right
        self.exact = exact
        if intcoords:
            if self.exact:
                self.gap = int(round(gap))
            else:
                self.gap = int(math.ceil(gap))
        else:
            self.gap = gap

    def rotateCW(self, n):
        """
        :param n: an integer in the set {1, 2, 3}
        :return: nothing

        We rotate the constraint 90*n degrees clockwise.
        """
        assert(1 <= n <= 3)
        dim = self.dim
        if (dim == adg.XDIM and n == 3) or (dim == adg.YDIM and n == 1):
            # These are the two cases in which the left and right shapes swap.
            self.left, self.right = self.right, self.left
        if n != 2:
            self.dim = 1 - dim

    def __repr__(self):
        s = ''
        d = {adg.XDIM: 'x', adg.YDIM: 'y'}[self.dim]
        s += 'SepCo: %s%s %s %s %s %s%s' % (
            self.left.ID, d,
            '+' if self.gap >= 0 else '-',
            abs(self.gap),
            '==' if self.exact else '<=',
            self.right.ID, d
        )
        return s

    def violation(self):
        """
        :return: the current violation of this separation constrain
        """
        lz = self.left.x if self.dim == adg.XDIM else self.left.y
        rz = self.right.x if self.dim == adg.XDIM else self.right.y
        dz = rz - lz
        vio = self.gap - dz
        if vio < 0:
            vio = -vio if self.exact else 0
        return vio

    def buildCCs(self, rs=None, ix=None):
        il = self.left.ID if ix is None else ix(self.left.ID)
        ir = self.right.ID if ix is None else ix(self.right.ID)
        sep = adg.SeparationConstraint(
            self.dim, il, ir, self.gap, self.exact
        )
        return [sep]

    def writeCpp(self, name='sep', ix=None):
        il = self.left.ID if ix is None else ix(self.left.ID)
        ir = self.right.ID if ix is None else ix(self.right.ID)
        s = ''
        s += '%s = new SeparationConstraint(%s, %d, %d, %.2f, %s);\n' % (name, {
            adg.XDIM: 'XDIM',
            adg.YDIM: 'YDIM'
        }[self.dim], il, ir, self.gap, str(self.exact).lower()
        )
        return s

    def writeDunnartSVG(self, nextID):
        # TODO: handle self.exact; set the position of the handles
        s = ''
        # Dunnart applies sepcos to guidelines, not shapes directly, so we
        # begin by creating an "alignment of one" for each of the two shapes
        # involved in the separation.
        for u in [self.left, self.right]:
            ac = AlignCo(self.dim, shapes=[[u, 0]])
            t, nextID = ac.writeDunnartSVG(nextID)
            s += t
        # Now just add the sepco svg on the two guidelines.
        direc = {adg.XDIM: 100, adg.YDIM: 101}[self.dim]
        relType = 'distribution' if self.exact else 'separation'
        s += '<dunnart:node type="%s" direction="%d" id="%d" ' % (
            relType, direc, nextID
        )
        s += 'sepDistance="%.2f" position="0">\n' % self.gap
        nextID += 1
        s += '    <dunnart:node type="constraint" relType="%s" ' % relType
        s += 'constraintID="%d" objOneID="%d" objTwoID="%d"/>\n' % (
            nextID - 1, nextID - 3, nextID - 2
        )
        s += '</dunnart:node>\n'
        return (s, nextID)

class DistCo:

    def __init__(self, dim, nodes, gap, flexible=False):
        self.dim = dim
        self.nodes = nodes
        self.gap = gap
        self.flexible = flexible

    def __repr__(self):
        s = ''
        s += 'DistCo %s %.2f:' % (
            {adg.XDIM: '--', adg.YDIM: '|'}[self.dim],
            self.gap
        )
        for node in self.nodes:
            s += ' %d' % node.ID
        return s

    def buildCCs(self, rs=None, ix=None):
        dist = adg.DistributionConstraint(self.dim)
        dist.isFlexible = self.flexible
        acs = []
        for u in self.nodes:
            ac = adg.AlignmentConstraint(self.dim)
            i = u.ID if ix is None else ix(u)
            ac.addShape(i, 0)
            acs.append(ac)
        for ac1, ac2 in zip(acs[:-1], acs[1:]):
            dist.addAlignmentPair(ac1, ac2)
        ccs = [dist] + acs
        return ccs

    def writeCpp(self, name='dist', ix=None):
        s = ''
        # TODO
        return s

    def rotateCW(self, n):
        """
        :param n: an integer in the set {1, 2, 3}
        :return: nothing

        We rotate the constraint 90*n degrees clockwise.
        """
        assert(1 <= n <= 3)
        dim = self.dim
        if (dim == adg.XDIM and n == 3) or (dim == adg.YDIM and n == 1):
            # These are the two cases in which the order of the shapes reverses.
            self.nodes.reverse()
        if n != 2:
            self.dim = 1 - dim


class FlexDistCo:

    def __init__(self, dim, p, q, r):
        """
        :param dim: adg.XDIM or adg.YDIM
        :param p: a Node
        :param q: a Node
        :param r: a Node
        """
        self.dim = dim
        self.p, self.q, self.r = p, q, r

    def rotateCW(self, n):
        assert(1 <= n <= 3)
        dim = self.dim
        if (dim == adg.XDIM and n == 3) or (dim == adg.YDIM and n == 1):
            # These are the two cases in which the order of the shapes reverses.
            self.r, self.p = self.p, self.r
        if n != 2:
            self.dim = 1 - dim

    def buildCCs(self, rs=None, ix=None):
        p = self.p if ix is None else ix(self.p)
        q = self.q if ix is None else ix(self.q)
        r = self.r if ix is None else ix(self.r)
        fdt = adg.FlexDistTripleConstraint(self.dim, p, q, r)
        return [fdt]

class FixedRelCo:

    def __init__(self, nodes):
        """
        :param nodes: a dict of Nodes by ID, to be constrained rigidly,
                      relative to one another
        """
        self.nodes = nodes

    def rotateCW(self, n):
        # Nothing to do.
        pass

    def buildCCs(self, rs=None, ix=None):
        if rs is None:
            raise Exception("FixedRelCo needs rs to build CC.")
        ixs = adg.Unsigneds()
        for node in self.nodes.values():
            ixs.push_back(ix(node))
        frc = adg.FixedRelativeConstraint(rs, ixs)
        return [frc]

class ProjSeq:

    def __init__(self):
        # sequence of sets of Python constraint objects:
        # These sets should be monotonic by dimension; i.e. each one should contain
        # the last one that was in the same dimension.
        self.PCsets = []
        # sequence of adg dimensions:
        self.dims = []
        # sequence of stress changes:
        self.dSes = []
        # pointer to next constraint set to be applied:
        self.ptr = 0
        # for maintaining monotonicity:
        self.finalSets = {adg.XDIM: [], adg.YDIM: []}

    def __repr__(self):
        s = ''
        s += 'ProjSeq:\n'
        for k in range(len(self.PCsets)):
            s += 'Set %d:\n' % k
            pcs = self.PCsets[k]
            for c in pcs:
                s += '    %s\n' % repr(c)
        return s

    def addConstraintSet(self, pcs, dim):
        """
        :param pcs: a list of SepCo and/or AlignCo objects
        :param dim: an adg dimension

        We add the new set of constraints, ensuring monotonicity by uniting with
        previous set of constraints in the same dimension, if any.
        """
        # Ensure monotonicity:
        finalSet = self.finalSets[dim]
        pcs = finalSet + pcs
        self.finalSets[dim] = pcs
        # Add the new set.
        self.PCsets.append(pcs)
        self.dims.append(dim)

    def __add__(self, other):
        ps = self.copy()
        qs = other.copy()
        qs.reset()
        for pcs, dim in qs:
            ps.addConstraintSet(pcs, dim)
        return ps

    def copy(self):
        ps = ProjSeq()
        ps.PCsets.extend(self.PCsets)
        ps.dims.extend(self.dims)
        ps.dSes.extend(self.dSes)
        ps.ptr = self.ptr
        ps.finalSets[adg.XDIM] = self.finalSets[adg.XDIM]
        ps.finalSets[adg.YDIM] = self.finalSets[adg.YDIM]
        return ps

    ######################################################################
    # Iteration
    #
    # We manage iteration over this object in such a way that it can be used
    # progressively, as new constraint sets are added to it.
    #
    # If, for example, the object has three constraint sets A, B, C, and you
    # iterate over these, then add two new sets D, E, and then iterate again,
    # the second iteration will be only over D and E, unless you call the
    # 'reset' method first.

    def __iter__(self):
        return self

    def next(self):
        z = zip(self.PCsets, self.dims)
        n = len(z)
        i = self.ptr
        if i >= n:
            raise StopIteration
        else:
            self.ptr = i + 1
            return z[i]

    def reset(self):
        self.ptr = 0

    ######################################################################

    def addStressChange(self, dS):
        self.dSes.append(dS)

    def getTotalStressChange(self):
        return sum(self.dSes)
