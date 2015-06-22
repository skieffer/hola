#!/usr/bin/python

import sys, allQAs

class Nbr:
    """
    Represents a neighbouring node to a central node.
    """

    def __init__(self, ID, dx, dy):
        """
        ID:    a unique identifier; useful for associating an
               instance with a particular node in the graph
               you are working on
        dx,dy: relative coordinates of this neighbour w.r.t.
               central node
        """
        self.ID = ID
        self.x = dx
        self.y = dy
        if dx == 0 and dy == 0:
            sys.stderr.write('ERROR: Nbr at (0,0).\n')

    def __repr__(self):
        return '%s'%self.ID

    def octalCode(self):
        """
        Semiaxes 0,1,2,3 get octal codes 0,2,4,6;
        Quadrants 0,1,2,3 get octal codes 1,3,5,7.
        """
        o = -1
        x = self.x; y = self.y
        if x > 0:
            o=7 if y<0 else 0 if y==0 else 1
        elif x == 0:
            o=6 if y<0 else 2
        else: # x < 0
            o=5 if y<0 else 4 if y==0 else 3
        return o

    def deflection(self):
        """
        If this nbr lies in quadrant n or on semiaxis n, then
        return the squared sine of the angle that this nbr makes
        with semiaxis n.
        Intuitively, this represents "how far into its quadrant
        this nbr lies."
        """
        x = self.x; y = self.y
        x2 = x*x; y2 = y*y
        l2 = float(x2 + y2)
        o = self.octalCode()
        defl = y2/l2 if o in [0,1,4,5] else x2/l2
        return defl


class Quad:
    """
    Represents a quadrant relative to a central node.
    """

    def __init__(self,num):
        """
        num: the number of this quadrant, 0, 1, 2, or 3
        """
        self.num = num
        self.nbrs = []
        self.costs = [0,0,0,0]

    def __repr__(self):
        s = ''
        s += 'Quad %d:\n'%self.num
        s += '    %s\n'%self.nbrs
        s += '    %s\n'%self.costs
        return s

    def addNbr(self,nbr):
        self.nbrs.append(nbr)

    def size(self):
        return len(self.nbrs)

    def sortAndComputeCosts(self):
        """
        To be called after all nbrs have been added.
        Sorts the nbrs into clockwise order, and then computes
        and stores the costs associated with each of the four
        possible actions D, A, B, and C.
        """
        if len(self.nbrs)==0:
            self.costs = [0,0,0,0]
            return
        self.nbrs = sorted(self.nbrs,
            key=lambda nbr: nbr.deflection()
        )
        A = self.nbrs[0].deflection()
        C = 1 - self.nbrs[-1].deflection()
        B = A + C
        self.costs = [0,A,B,C]

    def assignment(self,i):
        """
        i: an index 0, 1, 2, or 3 representing one of the
           four possible actions D, A, B, or C, respectively.
        We return an Assignment indicating which Nbr(s) were assigned
        to which semiaxis, and representing the cost of the
        assignment.
        """
        i = int(i)
        semis = [[], [], [], []]
        # j will be semiaxis to which nbr 0 should be assigned, if any;
        # k will be semiaxis to which nbr -1 should be assigned, if any.
        j = -1; k = -1
        if i in [1,2]: # The action is A or B
            j = self.num
        if i in [2,3]: # The action is B or C
            k = (self.num + 1) % 4
        if j >= 0:
            semis[j].append(self.nbrs[0])
        if k >= 0:
            semis[k].append(self.nbrs[-1])
        cost = self.costs[i]
        a = Assignment(semis, cost)
        return a

class Perm:
    """
    A permutation in D4.
    """

    def __init__(self,code):
        """
        code: a vector [a,b,c,d] indicating that this permutation
              should map 0 to a, 1 to b, 2 to c, and 3 to d.
        """
        self.lookup = code

    def __repr__(self):
        return str(self.lookup)

    def __call__(self,a):
        if type(a)==type(0):
            # If it is an integer, simply say where this permutation
            # maps that integer. (Must be in {0,1,2,3}.)
            return self.lookup[a]
        elif type(a)==type([]):
            # If it is a vector, act on the vector, i.e. permute it.
            b = [9,9,9,9]
            for j,k in zip(self.lookup,a):
                b[j] = k
            return b

    def inv(self):
        """
        Return a Perm object which represents the inverse
        permutation to this one.
        """
        code = [9,9,9,9]
        L = self.lookup
        for i in range(4):
            code[L[i]] = i
        I = Perm(code)
        return I

    def isFlip(self):
        """
        Say whether this permutation is one of the four flips in D4.
        This of course is the case just when the numbers occur in
        decreasing order, from left to right (mod 4).
        """
        a = self.lookup[0]; b = self.lookup[1]
        return (a-b)%4 == 1

def swapDir(qa):
    """
    If qa represents a quad action (i.e. is a 4-vector of
    integers 0,1,2,3), then "swap its direction", i.e. change
    all 1's to 3's, and vice versa.
    """
    qap = []
    for a in qa:
        if a%2==1: a = 4-a
        qap.append(a)
    return qap


class Assignment:
    """
    A struct to represent an assignment of nbrs to
    semiaxes, and the cost of this assignment (from some
    starting point).
    """

    def __init__(self, semis, cost):
        """
        semis: a list [t,u,v,w] of lists of Nbrs to be assigned to the
              semiaxes S0, S1, S2, S3, respectively.
              Entries may be empty lists, indicating that no nbr
              is to be assigned to that semiaxis.
        cost: the cost of this assignment
        """
        self.semis = semis
        self.cost = cost

    def __repr__(self):
        s = ''
        data = (self.card(),) + self.prettyIDTuple() + (self.cost,)
        s += 'Asgn %d:   %s %s %s %s   (%.3f)'%data
        return s

    def prettyIDTuple(self):
        return tuple([
            str([repr(nbr) for nbr in semiaxis])
            for semiaxis in self.semis
        ])

    def card(self):
        """
        'cardinality' i.e. number of assignments made
        """
        return sum([
            len(semiaxis) for semiaxis in self.semis
        ])

    def union(self, other):
        """
        Return a new Assignment representing the union
        of this one with the other.

        The semiaxis lists are merged, and the costs are summed.
        """
        return Assignment(
            [   list(set(s).union(set(o)))
                for s,o in zip(self.semis, other.semis)
            ],
            self.cost + other.cost
        )

def nullAssignment():
    a = Assignment([[], [], [], []], 0)
    return a


class Arrangement:
    """
    Represents the arrangement of nodes around a given centre
    node c.
    """

    def __init__(self, nbrs):
        """
        nbrs: list of Nbr objects representing the neighbouring
              nodes
        """

        #IDs = [nbr.ID for nbr in nbrs]
        #if 9 in IDs and 10 in IDs and 0 in IDs:
        #    print 'foo'

        self.nbrs = nbrs
        self.semis = [[], [], [], []]
        self.quads = [ Quad(0), Quad(1), Quad(2), Quad(3) ]
        # Assign nbrs to semis and quads based on their octals.
        for nbr in self.nbrs:
            o = nbr.octalCode()
            if o%2==0:
                s = o/2
                self.semis[s].append(nbr)
            else:
                q = (o-1)/2
                self.quads[q].addNbr(nbr)
        # Now all quads can sort and compute costs.
        for Q in self.quads:
            Q.sortAndComputeCosts()
            #print Q

    def vacancy(self):
        """
        Return a vector of nonnegative integers indicating how many
        neighbours are assigned to each semiaxis.
        """
        return [(0 if len(s)==0 else 1) for s in self.semis]

    def dist(self):
        """
        Compute and return the distribution vector.
        """
        return [Q.size() for Q in self.quads]

    def rDist(self):
        """
        Compute and return the reduced distribution.
        """
        d = self.dist()
        for i in range(4):
            if d[i] > 2: d[i] = 2
        return d

    def oriRedDistAndPhi(self):
        """
        Let d be the reduced distribution.
        Compute and return an ordered pair (d',phi) being
        the oriented reduced distribution d' and the orienting
        permutation phi satisfying d' = phi d.
        """
        # Get reduced distribution.
        d = self.rDist()
        # Initialise f.
        # By applying the same transformations to f as to d,
        # we will compute phi^-1 satisfying phi^-1 d' = d.
        f = [0,1,2,3]
        # We must also compute the permutation of the semiaxes.
        s = [0,1,2,3]
        # Compute minimum, and rotate so that p = m.
        m = min(d)
        i0 = d.index(m)
        d = d[i0:]+d[:i0]
        f = f[i0:]+f[:i0]
        s = s[i0:]+s[:i0]
        # If q > s, flip over main diagonal.
        if d[1] > d[3]:
            d = [d[j] for j in [0,3,2,1]]
            f = [f[j] for j in [0,3,2,1]]
            s = [s[k] for k in [1,0,3,2]]
        # If q = p and r > s, flip over vertical axis.
        if d[1] == d[0] and d[2] > d[3]:
            d = [d[j] for j in [1,0,3,2]]
            f = [f[j] for j in [1,0,3,2]]
            s = [s[k] for k in [2,1,0,3]]
        phi = Perm(f).inv()
        sigma = Perm(s).inv()
        return (d,phi,sigma)

    def basicAssignment(self):
        return Assignment(self.semis, 0)

    def getAssignmentForQuadAction(self, qa):
        """
        qa: a "quad action," i.e. a list [w,x,y,z], whose entires
            are integers 0, 1, 2, or 3, representing the quadrant
            actions D, A, B, C, respectively, with w, x, y, z
            intended for quadrants Q0, Q1, Q2, Q3 respectively.

        We return an Assignment object representing the semiaxis
        assignments and cost that correpond to the given quad action.
        """
        #print 'Getting assignment for QA: %s'%qa
        a = self.basicAssignment()
        for Q,i in zip( self.quads, qa ):
            ap = Q.assignment(i)
            a = a.union(ap)
        return a

    def getAsgns(self):
        """
        Return a list of Assignment objects indicating all
        possible 4- and 3-assignments around this centre node,
        sorted in order of preferability.

        The ordering is such that any and all 4-assignments
        come first, in order of ascending cost, followed by
        any and all 3-assignments, again in order of
        ascending cost.
        """
        trials = []
        trials.extend(self.listTrialNAssigns(4))
        trials.extend(self.listTrialNAssigns(3))
        trials.extend(self.listTrialNAssigns(2))
        trials.extend(self.listTrialNAssigns(1))
        return trials

    def listTrialNAssigns(self,N):
        """
        To have an N-assignment means that precisely N of the
        compass directions around the centre node have nbrs assigned
        to them.

        Return a list of Assignment objects indicating all
        possible N-assignments around this centre node,
        sorted in order of ascending cost.
        """
        # There are only four compass directions, so N cannot
        # exceed 4.
        if N > 4: return []
        # We need to have at least as many neighbours as the number
        # of assignments we are asked to make.
        if len(self.nbrs) < N: return []
        # Compute trial N-assignments.
        # Start with the vacancy vector.
        vac = self.vacancy()
        # Some semiaxes may already have Nbrs assigned to them.
        # Let n be number of semiaxes waiting to be filled
        # before we have an N-assignment.
        n = N - sum(vac)
        # And let f be the number of free neighbours, i.e. those in
        # the quads.
        f = sum([q.size() for q in self.quads])
        # If n == 0, then there is only one possible N-assignment,
        # namely, the one that is already there.
        if n == 0:
            a = self.basicAssignment()
            return [a]
        elif n < 0:
            # In this case we already have /more/ than N semiaxes
            # with neighbours on them, so an N-assignment,
            # strictly speaking, is not possible.
            return []
        elif f < n:
            # There are n > 0 semiaxes waiting to be filled before
            # we have an N-assignment, but this node has fewer than
            # n free nbrs, so there is no way to make an N-assignment.
            return []
        else:
            # There are n > 0 semiaxes waiting to be filled before
            # we have an N-assignment, and this node has at least
            # n free nbrs, so it is potentially possible to make
            # an N-assignment by assigning one or more further nbrs
            # to compass dirs.
            # Get the oriented reduced distribution
            # and the orienting permutation.
            d, phi, sigma = self.oriRedDistAndPhi()
            # Get the list of all possible quadrant actions.
            qas = listAllQuadActions(n,vac,d,phi,sigma)
            trials = []
            for qa in qas:
                a = self.getAssignmentForQuadAction(qa)
                trials.append(a)
            trials.sort(key=lambda a: a.cost)
            return trials

def listAllQuadActions(n,vac,d,phi,sigma):
    """
    n: the number of assignments to be made
    vac: a binary vector where 0 means a vacant semiaxis
    d: the oriented reduced distribution
    phi: the orienting permutation
    sigma: the corresponding permutation on the semiaxes

    We return a list of all possible quad actions which assign
    precisely n nodes to the vacant semiaxes, for the /original/
    distribution (which equals d iff phi is the identity perm).
    """
    # Write string description ds of d, and access allQAs library
    # for all n-assignments for d.
    ds = ''.join([str(k) for k in d])
    qasByCode = allQAs.qas[ds][n]
    #print 'For ds,n = %s,%s, get:\n%s\n'%(ds,n,qasByCode)
    # qasByCode is a dictionary whose keys are binary codes representing
    # which semiaxes are used, and whose values are lists of quad
    # action codes.
    # Apply sigma to the vacancy vector.
    vacp = sigma(vac)
    #print 'vacp: %s'%vacp
    # Make binary code for the permuted vacancy vector.
    s = 1
    vacpCode = 0
    for v in vacp:
        vacpCode += (1-v)*s
        s *= 2
    #print 'vacpCode: %s'%vacpCode
    # Build list of acceptable quad actions.
    q = []
    # We will need to permute them by phi inverse so that they
    # apply to the original distribution.
    phiInv = phi.inv()
    #print 'phi^-1 = %s'%phiInv
    for code in qasByCode:
        # qas[code] is acceptable iff the semiaxes that it uses
        # are a subset of those which are vacant.
        if code & vacpCode == code:
            orientedQAs = qasByCode[code]
            # Convert strings to lists.
            orientedQAs = [ [int(a) for a in qa] for qa in orientedQAs]
            #print 'orientedQAs: %s'%orientedQAs
            qas = map(phiInv, orientedQAs)
            # If phi^-1 is a flip, then must also swap 1's and 3's.
            if phiInv.isFlip():
                qas = map(swapDir, qas)
            #print 'qas: %s'%qas
            q.extend(qas)
    return q


def getAssignmentsForNode(node):
    """
    :param node: a Node object as defined in the pydaptagrams graphs module
    :return: all the possible assignments of the neighbours of this node to
             the compass directions. We return a list of Assignment objects,
             sorted in order of descending desirability.
    """
    nbrs = []
    for edge in node.edges.values():
        other = edge.otherEnd(node)
        dx = other.x - node.x
        dy = other.y - node.y
        nbr = Nbr(other.ID, dx, dy)
        nbrs.append(nbr)
    arr = Arrangement(nbrs)
    asgns = arr.getAsgns()
    return asgns

######################################################################
# Tests

def test1():
    """
    One in each of quads 0, 1, 2; two in quad 3.
    Should be a unique best assignment.
    """
    print 'Test 1'
    n0 = Nbr(0, 0.5, 0.866)
    n1 = Nbr(1, -1, 0.707)
    n2 = Nbr(2, -0.707, -0.707)
    n3 = Nbr(3, 0.5, -0.866)
    n4 = Nbr(4, 0.866, -0.5)
    nbrs = [n0, n1, n2, n3, n4]
    testNbrs(nbrs)

def test2():
    """
    Same as test1 except that now n2 is already assigned
    to semiaxis 2.
    """
    print 'Test 2'
    n0 = Nbr(0, 0.5, 0.866)
    n1 = Nbr(1, -1, 0.707)
    n2 = Nbr(2, -1, 0)
    n3 = Nbr(3, 0.5, -0.866)
    n4 = Nbr(4, 0.866, -0.5)
    nbrs = [n0, n1, n2, n3, n4]
    testNbrs(nbrs)

def test3():
    """
    All four are already assigned. There should be nothing
    else to do.
    """
    print 'Test 3'
    n0 = Nbr(0, 1, 0)
    n1 = Nbr(1, 0, 1)
    n2 = Nbr(2,-1, 0)
    n3 = Nbr(3, 0,-1)
    nbrs = [n0, n1, n2, n3]
    testNbrs(nbrs)

def test4():
    """
    Example demonstrating that a greedy approach wouldn't work.
    """
    print 'Test 4'
    n0 = Nbr(0, 0.5, 0.866)
    n1 = Nbr(1, -0.866, 0.5)
    n2 = Nbr(2,-1, -0.5)
    n3 = Nbr(3, 0.866,-0.5)
    nbrs = [n0, n1, n2, n3]
    testNbrs(nbrs)

def test5():
    """
    The original motivating example -- a case where the most
    preferable 3-assignment can fail due to an existing alignment
    constraint, and we need to choose the best Plan B.

           v
               c   w
        u          x

    Alignments: H(c,w); V(w,x); H(u,x)

    If a and b are the angles that uc and vc make with semiaxis 2,
    respectively, then assume we have 0 < a < b < 45.
    In this case the assignment (w - u v) is the most preferred,
    but it is thwarted by the existing alignments.
    Our algorithm successfully identifies (w u v -) as the next
    best alternative.
    """
    print 'Test 5'
    n0 = Nbr(0, 1, 0)
    n1 = Nbr(1, -1, 0.5)
    n2 = Nbr(2, -0.866, -0.5)
    nbrs = [n0, n1, n2]
    testNbrs(nbrs)

def testNbrs(nbrs):
    arr = Arrangement(nbrs)
    asgns = arr.getAsgns()
    for asgn in asgns:
        print repr(asgn)
    print

if __name__ == '__main__':
    test1()
    test2()
    test3()
    test4()
    test5()

