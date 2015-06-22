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

def text(x, y, t, style={}):
    s = ''
    if 'centre' in style:
        y += 5
        n = len(t)
        x -= {
            1: 5,
            2: 9,
            3: 14
        }[n]
    s += '<text x="%s" y="%s"' % (
        x, y
    )
    default = {
        'font-family': 'sans-serif'
    }
    for k, v in default.items():
        val = style.get(k, default[k])
        s += ' %s="%s"' % (k, val)
    s += '>%s</text>' % t
    return s

def rect(x, y, w, h, style={}):
    s = ''
    s += '<rect x="%s" y="%s" width="%s" height="%s"' % (
        x, y, w, h
    )
    default = {
        'stroke': 'black',
        'stroke-width': 1,
        'fill': 'none',
        'onclick': 'null',
        'id': ''
    }
    for k, v in default.items():
        val = style.get(k, v)
        s += ' %s="%s"' % (k, val)
    s += '/>'
    return s

def rectC(cx, cy, w, h, style={}):
    x, y = cx - w/2.0, cy - h/2.0
    return rect(x, y, w, h, style)

def rectBox(x, X, y, Y, style={}):
    w, h = X - x, Y - y
    return rect(x, y, w, h, style)


def line(x1, y1, x2, y2, style={}):
    s = ''
    s += '<line x1="%s" y1="%s" x2="%s" y2="%s"' % (
        x1, y1, x2, y2
    )
    default = {
        'stroke': 'black',
        'stroke-width': 1
    }
    for k, v in default.items():
        val = style.get(k, v)
        s += ' %s="%s"' % (k, val)
    s += '/>'
    return s

def path(d, style={}):
    s = ''
    s += '<path d="%s"' % d
    default = {
        'stroke': 'black',
        'stroke-width': 1,
        'fill': 'none'
    }
    for k, v in default.items():
        val = style.get(k, v)
        s += ' %s="%s"' % (k, val)
    s += '/>'
    return s

def roundedOrthoConnectorData(pts, curveRadius=10):
    d = ''
    if len(pts) < 2: return d
    cr0 = curveRadius
    p = pts[0]
    #d += 'M '+p[0]+','+p[1]
    d += 'M %.2f,%2f' % tuple(p)
    for i in range(1, len(pts) - 1):
        q = pts[i]
        r = pts[i+1]
        if q[1] < p[1]:
            #   q
            #   |
            #   p
            a = abs(p[1]-q[1])
            b = abs(q[0]-r[0])
            cr = min(cr0,a/2,b/2)
            #d += ' V '+(q[1]+cr)
            d += ' V %.2f' % (q[1]+cr)
            if r[0] < q[0]:
                # r -- q
                #      |
                #      p
                #d += ' a '+cr+','+cr+' 0 0,0 '+(-cr)+','+(-cr)
                d += ' a %.2f,%.2f 0 0,0 %.2f,%.2f' % (cr, cr, -cr, -cr)
            else:
                # q -- r
                # |
                # p
                #d += ' a '+cr+','+cr+' 0 0,1 '+(cr)+','+(-cr)
                d += ' a %.2f,%.2f 0 0,1 %.2f,%.2f' % (cr, cr, cr, -cr)
        elif q[1] > p[1]:
            #   p
            #   |
            #   q
            a = abs(p[1]-q[1])
            b = abs(q[0]-r[0])
            cr = min(cr0,a/2,b/2)
            #d += ' V '+(q[1]-cr)
            d += ' V %.2f' % (q[1]-cr)
            if (r[0] < q[0]):
                #      p
                #      |
                # r -- q
                #d += ' a '+cr+','+cr+' 0 0,1 '+(-cr)+','+(cr)
                d += ' a %.2f,%.2f 0 0,1 %.2f,%.2f' % (cr, cr, -cr, cr)
            else:
                # p
                # |
                # q -- r
                #d += ' a '+cr+','+cr+' 0 0,0 '+(cr)+','+(cr)
                d += ' a %.2f,%.2f 0 0,0 %.2f,%.2f' % (cr, cr, cr, cr)
        elif (q[0] > p[0]):
            #
            # p -- q
            #
            a = abs(p[0]-q[0])
            b = abs(q[1]-r[1])
            cr = min(cr0,a/2,b/2)
            #d += ' H '+(q[0]-cr)
            d += ' H %.2f' % (q[0]-cr)
            if (r[1] < q[1]):
                #      r
                #      |
                # p -- q
                #d += ' a '+cr+','+cr+' 0 0,0 '+(cr)+','+(-cr)
                d += ' a %.2f,%.2f 0 0,0 %.2f,%.2f' % (cr, cr, cr, -cr)
            else:
                # p -- q
                #      |
                #      r
                #d += ' a '+cr+','+cr+' 0 0,1 '+(cr)+','+(cr)
                d += ' a %.2f,%.2f 0 0,1 %.2f,%.2f' % (cr, cr, cr, cr)
        elif (q[0] < p[0]):
            #
            # q -- p
            #
            a = abs(p[0]-q[0])
            b = abs(q[1]-r[1])
            cr = min(cr0,a/2,b/2)
            #d += ' H '+(q[0]+cr)
            d += ' H %.2f' % (q[0]+cr)
            if (r[1] < q[1]):
                # r
                # |
                # q -- p
                #d += ' a '+cr+','+cr+' 0 0,1 '+(-cr)+','+(-cr)
                d += ' a %.2f,%.2f 0 0,1 %.2f,%.2f' % (cr, cr, -cr, -cr)
            else:
                # q -- p
                # |
                # r
                #d += ' a '+cr+','+cr+' 0 0,0 '+(-cr)+','+(cr)
                d += ' a %.2f,%.2f 0 0,0 %.2f,%.2f' % (cr, cr, -cr, cr)
        # Finally, set q as previous point, for next iteration.
        p = q
    p = pts[len(pts)-1]
    #d += 'L '+p[0]+','+p[1]
    d += 'L %.2f,%2f' % tuple(p)
    return d
