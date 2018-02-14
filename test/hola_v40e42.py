
import sys, os
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

import basic_tools
from hola.hola import HolaConfig
from hola.logging import LogLevel

def main():
    name = 'random/v40e42'
    config = HolaConfig()
    config.useFastSettings(True)
    config.LOG_LEVEL_GENERAL = LogLevel.TIMING
    config.LOG_LEVEL_TREE_PLACEMENT = LogLevel.TIMING
    #config.LOG_LEVEL_GENERAL = LogLevel.FINER_STAGE_GRAPHS
    #config.LOG_LEVEL_TREE_PLACEMENT = LogLevel.FINER_STAGE_GRAPHS
    basic_tools.doHOLA(name, config=config)

if __name__ == '__main__':
    main()

