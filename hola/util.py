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


class NearbyObjectFinder:
    """
    Suppose you are working with some objects A1, A2, ... each of which
    has a point (x, y) associated with it.
    For each new object Ai, you want to check whether there is already another
    object Aj whose coordinates are almost the same, within a given
    threshold.
    A NearbyObjectFinder can be used for this problem.
    Construct it with the desired threshold.
    Before adding any new object to it, use its 'findObject' method to
    see whether it already has an object with both x and y within the
    threshold of the x and y for the new object.
    If not, then add the object using the 'addObject' method.
    """

    def __init__(self, threshold):
        self.th = threshold
        self.objects = {}

    def addObject(self, x, y, obj):
        """
        Add a new object, and say what its x,y-coords are
        :param x: the x-coord for the object (float)
        :param y: the y-coord for the object (float)
        :param obj: the object
        :return: nothing
        """
        # We use buckets by storing the object under the rounded integer coordinates...
        xr, yr = int(round(x)), int(round(y))
        x_dict = self.objects.get(xr, {})
        # ...but we store it along with the given float coordinates.
        y_list = x_dict.get(yr, [])
        y_list.append((x, y, obj))
        x_dict[yr] = y_list

        self.objects[xr] = x_dict

    def findObject(self, x, y):
        """
        Check to see if any object has been stored yet with coordinates that are
        within the threshold of the given ones. If so, return it, else None.
        :param x: target x-coord (float)
        :param y: target y-coord (float)
        :return: object stored under nearby coords, else None
        """
        # Consider the range of all integer coords under which a nearby
        # point could have been stored, according to the threshold.
        th = self.th
        u0, u1 = int(round(x - th)), int(round(x + th))
        v0, v1 = int(round(y - th)), int(round(y + th))
        U = range(u0, u1+1)
        V = range(v0, v1+1)
        # Prepare return value.
        obj = None
        # Check to see if any object was stored with coordinates that are
        # within the threshold of the given ones.
        for u in U:
            if not self.objects.has_key(u): continue
            x_dict = self.objects[u]
            for v in V:
                if not x_dict.has_key(v): continue
                y_list = x_dict[v]
                for x1, y1, obj1 in y_list:
                    if abs(x - x1) < th and abs(y - y1) < th:
                        obj = obj1
                        return obj
        return obj
