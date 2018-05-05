#!/usr/bin/python

"""
Generate a module (write to stdout) encoding ALL the posible
assignments ("quad actions") for oriented reduced distributions.

A quad action is a sequence of four numerical codes, saying what to do in each
of the four quadrants, Q0, Q1, Q2, Q3, respectively. Before we switched to
numerical codes we used to use letter codes, and the mapping and meanings are
as follows:

    Num Char Meaning
    --- ---- -------
    0   D    "Do Nothing"
    1   A    "Anticlockwise": assign the node nearest the anticlockwise end
                              of the interval to that semiaxis. E.g. in Q0 the
                              node nearest semiaxis 0 gets assigned to it.
    2   B    "Both": do both the Clockwise and Anticlockwise assignments.
    3   C    "Clockwise": assign the node nearest the clockwise end of the
                          interval to that semiaxis. E.g. in Q0 the node
                          nearest semiaxis 1 gets assigned to it.

The output will look like this:

qas = {
    '0112':{    # oriented reduced distribution
        '4':{   # number of edges wanting an assignment
            '15':[    # which semiaxes are used, binary encoded
                '0112'  # quad actions
            ]
        },
        '3':{
            '7':[ # use semiaxes 012
                '0113'
            ],
            '11':[ # use semiazes 013
                '0102',
                '0133'
            ]
        },
        '2':{

        },
        '1':{

        }
    }
}
"""

import sys, json

# All 20 oriented reduced distributions:
ords = [
    (0,0,0,1),
    (0,0,0,2),
    (0,0,1,1),
    (0,1,0,1),
    (0,0,1,2),
    (0,1,0,2),
    (0,0,2,2),
    (0,2,0,2),
    (0,1,1,1),
    (0,1,1,2),
    (0,1,2,1),
    (0,1,2,2),
    (0,2,1,2),
    (0,2,2,2),
    (1,1,1,1),
    (1,1,1,2),
    (1,1,2,2),
    (1,2,1,2),
    (1,2,2,2),
    (2,2,2,2)
]

class Vect:

    def __init__(self,*args):
        self.comps = args

    def __add__(self,other):
        # We assume they are the same length.
        s = []
        for a,b in zip(self.comps,other.comps):
            s.append(a+b)
        w = Vect(0)
        w.comps = s
        return w

    def __repr__(self):
        s = ''
        for c in self.comps:
            s += ',%s'%c
        s = '('+s[1:]+')'
        return s

    def invalid(self):
        ans = False
        for c in self.comps:
            if c > 1:
                ans = True
                break
        return ans

    def height(self):
        return sum(self.comps)

    def binaryUsageCode(self):
        """
        Return an integer from 0 to 15 indicating which semiaxes
        are occupied, encoding semiaxis Sn in the nth binary bit.
        """
        code = 0
        s = 1
        for c in self.comps:
            if c > 0: code += s
            s *= 2
        return code

    def missing(self):
        """
        Return a list of the indices that get a 0, i.e. to which
        no node has been assigned.
        """
        m = []
        for i in range(len(self.comps)):
            if self.comps[i] == 0:
                m.append(i)
        return m

#    |   
#  r | s 
#    |
# ---+---
#    |
#  q | p
#    |

pH = [
    Vect(0,0,0,0),
    Vect(1,0,0,0),
    Vect(1,1,0,0),
    Vect(0,1,0,0)
]

qH = [
    Vect(0,0,0,0),
    Vect(0,1,0,0),
    Vect(0,1,1,0),
    Vect(0,0,1,0)
]

rH = [
    Vect(0,0,0,0),
    Vect(0,0,1,0),
    Vect(0,0,1,1),
    Vect(0,0,0,1)
]

sH = [
    Vect(0,0,0,0),
    Vect(0,0,0,1),
    Vect(1,0,0,1),
    Vect(1,0,0,0)
]


def enumerate(d):
    """
    d: an oriented reduced distribution

    We write the dictionary for this, enumerating all legal
    quad actions, as keyed by number of edges to be assigned
    (primarily) and vacant semiaxes (secondarily).
    """
    out = sys.stdout
    choices = []
    N = 1
    for c in d:
        if c==0:
            choices.append([0])
        elif c==1:
            choices.append([0,1,3])
            N *= 3
        elif c==2:
            choices.append([0,1,2,3])
            N *= 4
    qas = {
        1:{}, 2:{}, 3:{}, 4:{}
    }
    for p in choices[0]:
        h0 = Vect(0,0,0,0) # E, S, W, N
        h1 = h0 + pH[p]
        for q in choices[1]:
            h2 = h1 + qH[q]
            for r in choices[2]:
                h3 = h2 + rH[r]
                for s in choices[3]:
                    hits = h3 + sH[s]
                    if hits.invalid(): continue
                    H = hits.height()
                    if H == 0: continue
                    assert(1 <= H and H <= 4)
                    B = hits.binaryUsageCode()
                    L = qas[H].get(B,[])
                    qa = '%d%d%d%d'%(p,q,r,s)
                    try:
                        i = L.index(qa)
                    except:
                        L.append(qa)
                    qas[H][B] = L
    return qas
    

def test1():
    d = (1,1,1,2)
    e = enumerate(d)
    print d
    print json.dumps(e,indent=4)

def writeFullBlockIndented():
    block = {}
    for d in ords:
        e = enumerate(d)
        block[''.join([str(k) for k in d])] = e
    j = json.dumps(block,indent=4)
    sys.stdout.write('qas='+j)

def writeFullBlockCompact():
    block = {}
    for d in ords:
        e = enumerate(d)
        block[''.join([str(k) for k in d])] = e
    r = repr(block)
    sys.stdout.write('qas='+r)

if __name__=='__main__':
    #writeFullBlockIndented()
    writeFullBlockCompact()

