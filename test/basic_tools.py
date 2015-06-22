
import sys, os

sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

from hola import hola
import hola.logging as logging
from hola.gmlparse import buildGraph

def doHOLA(name, config=None):
    with open('graphs/%s.gml' % name) as f:
        gml = f.read()
    G = buildGraph(gml)
    groupFolder = 'output'
    parts = name.split('/')
    uniqueName = parts[-1]
    if len(parts) == 2: groupFolder += '/' + parts[0]
    final = groupFolder + '/final'
    if not os.path.exists(final):
        os.makedirs(final)
    logger = logging.Logger('output/'+name, uniqueName)
    projLogger = logging.Logger('output/'+name+'/proj', uniqueName)
    if config is None:
        config = hola.HolaConfig()
    # Do layout
    print '='*70
    print 'HOLA Layout for %s' % name
    hola.hola(G, config=config, logger=logger, projLogger=projLogger)
    cc = config.writeConfigCode()
    logger.writeGML("_final_layout_%s" % cc, graph=G)
    if logger.level >= logging.LogLevel.TIMING:
        td = logger.rootTimer.getTimingDataRec()
        print
        print 'TOTAL TIMES:'
        print td.write()
    if projLogger.level >= logging.LogLevel.TIMING:
        td = projLogger.rootTimer.getTimingDataRec()
        print
        print 'TOTAL TIMES FOR TREE PLACEMENT:'
        print td.write()
    # Copy final result
    os.system('cp output/%s/%s_final_layout*.gml %s' % (
        name, uniqueName, final
    ))

