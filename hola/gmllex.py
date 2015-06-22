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

"""
Lex a GML file and create a list of tokens.
"""

import sys, re


class Token:
    def __init__(self, type, attr=None):
        self.type = type
        self.attr = attr
        self.row = None
        self.col = None
    def __cmp__(self, o):
        return cmp(self.type, o)
    def __str__(self):
        s = 'Token: type=%s'%self.type
        if self.attr: s += ' attr=%s'%self.attr
        if self.row: s += ' row=%s'%self.row
        if self.col: s += ' col=%s'%self.col
        return s
    def __repr__(self):
        return str(self)
    def setRowCol(self, r, c):
        self.row = r
        self.col = c

_keywords = r"""
    \[
    \]
""".split()

# Longest match will be taken.
# If two or more token types match at the same length, then the one
# listed earliest in this list will be taken.
tokenTypes = [
    ('_keyword'  , re.compile('|'.join(_keywords)) ),
    ('key'       , re.compile(r'[A-Za-z]+') ),
    ('int'       , re.compile(r'-?\d+') ),
    ('float'     , re.compile(r'-?\d+\.\d+') ),
    ('string'    , re.compile(r'"[^"]*"') ),
    ('_skip'     , re.compile(r'\s+') )
]

def lex(text):

    class Stream:
        "Keeps track of a source text, row, and column."
        def __init__(self,text):
            self.text = text; self.row = 1; self.col = 1;
        def __len__(self):
            return len(self.text)
        def __getitem__(self,i):
            return self.text[i]
        def consumeChars(self,n):
            block = self.text[:n]
            self.text = self.text[n:]
            lines = block.split('\n')
            nls = len(lines)-1
            self.row += nls
            if nls == 0:
                self.col += len(block)
            else:
                self.col = 1 + len(lines[-1])
            return block
        def raiseUnclosedBlockException(self,tokens=None):
            err = ''
            if tokens:
                err += 'Tokens so far matched:\n'+repr(tokens)+'\n\n'
            err += 'Lexing error. Block beginning at '
            err += 'row %s, column %s '%(self.row,self.col)
            if len(self) > 0:
                n = min(len(self),50)
                err += '\n%s'%self.text[:n]
                if n < len(self): err += '...'
                err += '\n'
            err += 'never finished.\n'
            raise Exception(err)
        def raiseNoMatchException(self,tokens=None):
            err = ''
            if tokens:
                err += 'Tokens so far matched:\n'+repr(tokens)+'\n\n'
            err += 'Lexing error. Characters at '
            err += 'row %s, column %s '%(self.row,self.col)
            if len(self) > 0:
                n = min(len(self),50)
                err += '\n%s'%self.text[:n]
                if n < len(self): err += '...'
                err += '\n'
            err += 'do not match any token type.\n'
            raise Exception(err)

    def isSkipType(kind):
        return kind[:5] == '_skip'

    def tryRE(regexp,kind,stream,tokens):
        """
        Try to match regexp (of kind kind) at head of stream.
        If match, then append new token to tokens, update
        stream, and return True; else return False.
        """
        M = regexp.match(stream.text)
        if M:
            block = M.group()
            if not isSkipType(kind):
                if kind == '_keyword':
                    token = Token(block)
                else:
                    token = Token(kind,attr=block)
                tokens.append(token)
            stream.consumeChars(len(block))
            return True
        else:
            return False

    def matchLongestRE(stream,tokens):
        """
        Match the head of the stream against all regexps,
        accept the longest match, and return True, unless
        there were no matches at all.
        """
        kind = None; block = ''
        for k,r in tokenTypes:
            m = r.match(stream.text)
            if m:
                g = m.group()
                #sys.stderr.write('%s: %s'%(k,g))
                if len(g) > len(block):
                    kind = k
                    block = g
        if kind:
            if not isSkipType(kind):
                if kind == '_keyword':
                    token = Token(block)
                else:
                    token = Token(kind,attr=block)
                token.setRowCol(stream.row, stream.col)
                tokens.append(token)
            stream.consumeChars(len(block))
            return True
        else:
            return False

    tokens = []
    stream = Stream(text)
    while len(stream) > 0:
        # Try to match all the token types, and accept
        # the longest match, or return false if none matches.
        if matchLongestRE(stream,tokens):
            pass
        else:
            # Nothing matches!
            stream.raiseNoMatchException(tokens=tokens)
    return tokens


def test1():
    text = sys.stdin.read()
    tokens = lex(text)
    for t in tokens:
        sys.stdout.write('%s\n'%t)

if __name__ == '__main__':
    test1()

