#!/usr/bin/python

"""
Simple application of FD layout
"""

import sys, os
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

from hola.gmlparse import buildGraph

def FD(path, out_path):
    #out_path = path[:-4] + '_FD.gml'
    with open(path) as f:
        gml = f.read()
    G = buildGraph(gml)
    avgdim = G.computeAvgNodeDim()
    iel = 2 * avgdim
    op = False
    G.fdlayout(iel, op)
    output_type = out_path[-3:]
    with open(out_path, 'w') as f:
        if output_type == 'gml':
            f.write(G.writeGML())
        elif output_type == 'svg':
            f.write(G.writeDunnartSVG())

USAGE = """\
Usage:

    $ ./test_fd INPUTFILE OUTPUTFILE

Input file must have 'gml' extension.
Output file must have 'gml' or 'svg' extension.
"""

if __name__ == "__main__":
    try:
        path = sys.argv[1]
        assert os.path.exists(path) and path[-3:] == 'gml'
        out_path = sys.argv[2]
        assert out_path[-3:] in ['gml', 'svg']
    except:
        print USAGE
        sys.exit(1)
    FD(path, out_path)
