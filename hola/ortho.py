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

class Compass:
    EAST = 0
    SOUTH = 1
    WEST = 2
    NORTH = 3
    SE = 4
    SW = 5
    NW = 6
    NE = 7

    cwCards = [0, 1, 2, 3]
    acwCards = [0, 3, 2, 1]

    cwOrds = [4, 5, 6, 7]
    acwOrds = [4, 7, 6, 5]

    cwRose = [0, 4, 1, 5, 2, 6, 3, 7]
    acwRose = [0, 7, 3, 6, 2, 5, 1, 4]

    horizontal = [0, 2]
    vertical = [1, 3]

    increasing = [0, 1]
    decreasing = [2, 3]

    abbrev = {
        EAST: "E", SOUTH: "S", WEST: "W", NORTH: "N",
        SE: "SE", SW: "SW", NW: "NW", NE: "NE"
    }

    # Directions w.r.t which we must /increase/ the const
    # coord in order to move to the right:
    rightSidePlus = [NORTH, EAST]
    # Directions w.r.t which we must /decrease/ the const
    # coord in order to move to the right:
    rightSideMinus = [SOUTH, WEST]

    variableDimension = {
        EAST: adg.XDIM,
        SOUTH: adg.YDIM,
        WEST: adg.XDIM,
        NORTH: adg.YDIM
    }

    constantDimension = {
        EAST: adg.YDIM,
        SOUTH: adg.XDIM,
        WEST: adg.YDIM,
        NORTH: adg.XDIM
    }

    components = {
        SE: [SOUTH, EAST],
        SW: [SOUTH, WEST],
        NW: [NORTH, WEST],
        NE: [NORTH, EAST],
        # For convenience:
        EAST: [EAST],
        SOUTH: [SOUTH],
        WEST: [WEST],
        NORTH: [NORTH]
    }

    signs = {
        EAST: (1, 0),
        SE: (1, 1),
        SOUTH: (0, 1),
        SW: (-1, 1),
        WEST: (-1, 0),
        NW: (-1, -1),
        NORTH: (0, -1),
        NE: (1, -1)
    }

    libavoidVisibility = {
        EAST: adg.ConnDirRight,
        SOUTH: adg.ConnDirDown,
        WEST: adg.ConnDirLeft,
        NORTH: adg.ConnDirUp
    }

    @classmethod
    def cwClosedInterval(cls, d0, d1):
        """
        :param d0: a direction
        :param d1: another direction
        :return: list of all compass directions from d0 to d1 inclusive, going clockwise
        """
        rr = cls.cwRose + cls.cwRose
        i0 = rr.index(d0)
        i1 = i0 + rr[i0:].index(d1)
        return rr[i0:i1+1]

    @classmethod
    def acwClosedInterval(cls, d0, d1):
        """
        :param d0: a direction
        :param d1: another direction
        :return: list of all compass directions from d0 to d1 inclusive, going anticlockwise
        """
        return list(reversed(cls.cwClosedInterval(d1, d0)))

    @classmethod
    def cwRoseDistance(cls, d0, d1):
        """
        :param d0: a direction
        :param d1: another direction
        :return: the number of steps on the compass rose going clockwise from d0 to d1
        """
        return len(cls.cwClosedInterval(d0, d1)) - 1

    @classmethod
    def shortestRoseDistance(cls, d0, d1):
        """
        :param d0: a direction
        :param d1: another direction
        :return: the minimum number of steps on the compass rose from d0 to d1, going in either direction
        """
        return min(cls.cwRoseDistance(d0, d1), cls.acwRoseDistance(d0, d1))

    @classmethod
    def acwRoseDistance(cls, d0, d1):
        """
        :param d0: a direction
        :param d1: another direction
        :return: the number of steps on the compass rose going anticlockwise from d0 to d1
        """
        return len(cls.acwClosedInterval(d0, d1)) - 1

    @classmethod
    def sameDimension(cls, d0, d1):
        """
        :param d0: a cardinal Compass direction
        :param d1: a cardinal Compass direction
        :return: boolean saying if these directions are in the same dimension
        """
        return (d0 % 2) == (d1 % 2)

    @classmethod
    def perpendicular(cls, d0, d1):
        """
        :param d0: a cardinal Compass direction
        :param d1: a cardinal Compass direction
        :return: boolean saying if these directions are perpendicular to one another
        """
        return not cls.sameDimension(d0, d1)

    @classmethod
    def cardinalDirection(cls, p1, p2):
        """
        :param p1: either a Node object, or the coords (x1, y1) of a point
        :param p2: either a Node object, or the coords (x2, y2) of a point
        :return: the predominant cardinal direction from p1 to p2
        """
        try:
            x1, y1 = p1.x, p1.y
        except:
            x1, y1 = p1
        try:
            x2, y2 = p2.x, p2.y
        except:
            x2, y2 = p2
        dx, dy = x2 - x1, y2 - y1
        if abs(dy) <= abs(dx):
            return cls.EAST if x1 < x2 else cls.WEST
        else:
            return cls.SOUTH if y1 < y2 else cls.NORTH

    @classmethod
    def possibleCardinalDirections(cls, node1, node2):
        """
        :param node1: a Node
        :param node2: a Node
        :return: a list of the possible cardinal directions from node1 to node2,
                 if they were to be aligned non-aggressively
        """
        x1, y1 = node1.x, node1.y
        x2, y2 = node2.x, node2.y
        dx, dy = x2 - x1, y2 - y1
        dirs = []
        if dx > 0:
            dirs.append(cls.EAST)
        if dx < 0:
            dirs.append(cls.WEST)
        if dy > 0:
            dirs.append(cls.SOUTH)
        if dy < 0:
            dirs.append(cls.NORTH)
        return dirs

    @classmethod
    def getRotationFunction(cls, fromDir, toDir):
        # For now we only handle cardinal directions.
        if fromDir not in cls.cwCards or toDir not in cls.cwCards:
            raise Exception("only cardinal directions are currently handled")
        a, b = fromDir, toDir
        d = (b - a) % 4
        return [
            lambda v: (v[0], v[1]),
            lambda v: (-v[1], v[0]),
            lambda v: (-v[0], -v[1]),
            lambda v: (v[1], -v[0])
        ][d]

    @classmethod
    def flip(cls, direc):
        i0 = cls.cwRose.index(direc)
        return cls.cwRose[(i0+4)%8]

    @classmethod
    def cw90(cls, direc):
        i0 = cls.cwRose.index(direc)
        return cls.cwRose[(i0+2)%8]

    @classmethod
    def acw90(cls, direc):
        i0 = cls.cwRose.index(direc)
        return cls.cwRose[(i0-2)%8]

    @classmethod
    def rotateCW(cls, n, direc):
        i0 = cls.cwRose.index(direc)
        return cls.cwRose[(i0+n)%8]

    @classmethod
    def vectorSigns(cls, direc):
        """
        :param direc: a Compass direction
        :return: (xs, ys) where xs in {-1, 0, 1} represents the sign of
        the x-coordinate of a vector lying in the "octant" represented
        by direc, and likewise for ys. Here an "octant" is a semiaxis for
        a cardinal direction and an open quadrant for an ordinal direction.
        """
        return cls.signs[direc]

    @classmethod
    def vector(cls, direc, mag=1):
        """
        :param direc: a Compass direction, cardinal or ordinal
        :param mag: a float
        :return: a vector in the form [x, y] having the given magnitude and
                 pointing in the given direction
        """
        hsqrt2 = 0.7071067811865476
        v = {
            cls.EAST: [1, 0],
            cls.SOUTH: [0, 1],
            cls.WEST: [-1, 0],
            cls.NORTH: [0, -1],
            cls.SE: [hsqrt2, hsqrt2],
            cls.SW: [-hsqrt2, hsqrt2],
            cls.NW: [-hsqrt2, -hsqrt2],
            cls.NE: [hsqrt2, -hsqrt2]
        }[direc]
        return [mag*c for c in v]

class LineSegment:

    EPSILON = 0.1

    def __init__(self, p0, p1):
        """
        :param p0: a point (x0, y0)
        :param p1: a point (x1, y1)
        """
        self.p0 = p0
        self.p1 = p1
        x0, y0 = p0
        x1, y1 = p1
        #assert(x0 == x1 or y0 == y1)
        self.direc = Compass.cardinalDirection(p0, p1)
        self.varDim = Compass.variableDimension[self.direc]
        self.constDim = Compass.constantDimension[self.direc]
        if self.direc in Compass.vertical:
            self.z, self.w0, self.w1 = x0, y0, y1
        else:
            self.z, self.w0, self.w1 = y0, x0, x1
        self.wl, self.wh = min(self.w0, self.w1), max(self.w0, self.w1)
        self.length = self.wh - self.wl

    def __repr__(self):
        s = ''
        s += 'Seg: %s from %f to %f at %f' % (
            '|' if self.direc in Compass.vertical else '--',
            self.wl, self.wh, self.z
        )
        return s

    def constCoord(self):
        return self.z

    def lowCoord(self):
        return self.wl

    def highCoord(self):
        return self.wh

    def closedIntervalIncludesPt(self, pt):
        w = pt[0] if self.direc in Compass.horizontal else pt[1]
        return self.closedIntervalIncludesCoord(w)

    def openIntervalIncludesPt(self, pt):
        w = pt[0] if self.direc in Compass.horizontal else pt[1]
        return self.openIntervalIncludesCoord(w)

    def closedIntervalIncludesCoord(self, w):
        return self.wl <= w and w <= self.wh

    def openIntervalIncludesCoord(self, w):
        eps = 0.0001
        return self.wh - w > eps and w - self.wl > eps

    def closedIntervalIntersects(self, I):
        "Use this method if BOTH intervals are closed."
        a, b = I
        return b >= self.wl and self.wh >= a

    def openIntervalIntersects(self, J):
        "Use this method if EITHER interval is open."
        DEBUG = False
        a, b = J
        ans = b - self.wl > LineSegment.EPSILON and self.wh - a > LineSegment.EPSILON
        if DEBUG:
            print '  --OPEN interval check--'
            print '    b = %.100f' % b
            print '   wl = %.100f' % self.wl
            print '   wh = %.100f' % self.wh
            print '    a = %.100f' % a
        return ans

    def closedIntervalIntersection(self, I):
        if not self.closedIntervalIntersects(I):
            return None
        else:
            a, b = I
            return [max(a, self.wl), min(b, self.wh)]

    def openIntervalIntersection(self, J):
        if not self.openIntervalIntersects(J):
            return None
        else:
            a, b = J
            return [max(a, self.wl), min(b, self.wh)]

    def ptOnWhichSide(self, pt):
        z = pt[0] if self.direc in Compass.vertical else pt[1]
        return self.coordOnWhichSide(z)

    def coordOnWhichSide(self, z):
        z0 = self.z
        return -1 if z < z0 else (1 if z > z0 else 0)


def getSideOfBox(box, side):
    """
    :param box: (bounding) box in the form (x, X, y, Y) giving extreme coords
    :param side: cardinal Compass direction indicating one of the four sides of the box
    :return: LineSegment representing the desired side of the box
    """
    x, X, y, Y = box
    if side == Compass.EAST:
        p0 = (X, y)
        p1 = (X, Y)
    elif side == Compass.SOUTH:
        p0 = (x, Y)
        p1 = (X, Y)
    elif side == Compass.WEST:
        p0 = (x, y)
        p1 = (x, Y)
    else:
        assert side == Compass.NORTH
        p0 = (x, y)
        p1 = (X, y)
    return LineSegment(p0, p1)
