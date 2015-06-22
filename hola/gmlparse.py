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

import sys, os
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

from thirdparty.SPARK.spark import GenericASTBuilder, GenericASTTraversal
from thirdparty.SPARK.ast import AST as BasicAST
from gmllex import lex
from graphs import *

class AST(BasicAST):
    def __init__(self, type):
        BasicAST.__init__(self, type)
        self.attr = None
    def __repr__(self):
        return self.paddedRep('')
    def paddedRep(self,pad):
        s = ''
        s += pad+self.type
        if self.attr: s += ': %s'%self.attr
        s +='\n'
        for k in self._kids:
            s += k.paddedRep(pad+'  ')
        return s

######################################################################

class AutoGmlParser(GenericASTBuilder):

    def __init__(self, AST, start='itemList'):
        GenericASTBuilder.__init__(self, AST, start)

    def p_module(self, args):
        '''
            itemList ::= itemList item
            itemList ::= item
            item ::= key value
            value ::= string
            value ::= int
            value ::= float
            value ::= [ itemList ]
            value ::= [ ]
        '''
        pass

    def terminal(self, token):
        #
        #  Homogeneous AST.
        #
        rv = AST(token.type)
        rv.attr = token.attr
        return rv

    def nonterminal(self, type, args):

        # Flatten "list types".
        if type==args[0].type:
            args[0]._kids.append(args[1])
            return args[0]
        # Pass through "alternation types".
        # These are precisely those nonterminals which get
        # exactly 1 argument AND are NOT list types.
        if len(args) == 1 and not type[-4:]=='List':
            return args[0]

        return GenericASTBuilder.nonterminal(self, type, args)

def parse(tokens):
    parser = AutoGmlParser(AST)
    try:
        ast = parser.parse(tokens)
    except Exception as e:
        msg = 'Problem parsing gml file:\n    %s\n'%e
        raise Exception(msg)
    return ast

######################################################################
# Graph builder

class GraphBuilder(GenericASTTraversal):

    def __init__(self, ast):
        GenericASTTraversal.__init__(self, ast)
        self.graph = None
        # Begin tree traversal:
        self.postorder()

    def n_itemList(self, astn):
        astn.dct = {}
        for item in astn:
            # If the item is one of the types that gets repeated
            # then add it to a list in the AST node's dictionary astn.dct.
            reptypes = 'node edge point'.split()
            if item.key in reptypes:
                k = item.key
                # Get the name of the corresponding list.
                listname = {
                    'node': 'nodeList',
                    'edge': 'edgeList',
                    'point': 'Line'
                }[k]
                L = astn.dct.get(listname, [])
                L.append( item.val )
                astn.dct[listname] = L
            # Otherwise assume it is a unique type, and just add
            # the key-value pair directly to the AST node's dict.
            else:
                astn.dct[item.key] = item.val

    def n_item(self, astn):
        key = astn[0].attr
        astn.key = key
        if len(astn[1]) == 3:
            # In this case it is
            # value ::= [ itemList ]
            # and we take the dictionary of the item list as the value.
            val = astn[1][1].dct
            # For special cases of objects we want to build, construct them on
            # the dictionary in astn.val, and then replace astn.val with the
            # constructed object.
            if key == 'graph':
                G = Graph()
                G.initFromGML(val)
                astn.val = G
                # We expect only one graph (is that right?) so set it globally too:
                self.graph = G
            elif key == 'node':
                N = Node()
                N.initFromGML(val)
                astn.val = N
            elif key == 'edge':
                E = Edge()
                E.initFromGML(val)
                astn.val = E
            else:
                # If no special case then just stash the dictionary as
                # the value of this AST node.
                astn.val = val
        elif len(astn[1]) == 2:
            astn.val = None
        else:
            # Otherwise it is a constant (string, int, or float)
            a = astn[1].attr
            # Interpret value if possible.
            if astn[1].type == 'int':
                a = int(a)
            elif astn[1].type == 'float':
                a = float(a)
            elif astn[1].type == 'string':
                a = a[1:-1] # (chop off quotation marks)
            astn.val = a


######################################################################
# Test

def test1():
    text = sys.stdin.read()
    tokens = lex(text)
    ast = parse(tokens)
    stdout.write('%s\n'%ast)

def test2():
    text = sys.stdin.read()
    G = buildGraph(text)
    stdout.write('%s\n'%repr(G))


######################################################################
# API

def buildGraph(gml):
    """
    gml: the text from a .gml file.

    We build and return a Graph object.
    """
    tokens = lex(gml)
    ast = parse(tokens)
    gb = GraphBuilder(ast)
    G = gb.graph
    return G

##########################
if __name__ == '__main__':
    test2()

