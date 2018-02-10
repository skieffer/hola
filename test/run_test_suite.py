
import sys, os
sys.path.append(os.path.join(os.getcwd(), os.path.pardir))

import basic_tools
from hola.hola import HolaConfig
from hola.logging import LogLevel

def run_tests(test_suite):
    # Set up the configuration
    config = HolaConfig()
    config.useFastSettings(True)
    config.LOG_LEVEL_GENERAL = LogLevel.TIMING
    config.LOG_LEVEL_TREE_PLACEMENT = LogLevel.TIMING
    #config.LOG_LEVEL_GENERAL = LogLevel.FINER_STAGE_GRAPHS
    #config.LOG_LEVEL_TREE_PLACEMENT = LogLevel.FINER_STAGE_GRAPHS

    # Run the tests
    listing_path = os.path.join('graphs', test_suite)
    init_path = os.path.split(test_suite)[0]
    with open(listing_path) as f:
        lines = f.readlines()
    for line in lines:
        name = line.strip()
        # Skip empty lines
        if not name: continue
        # Skip lines beginning with "#"
        if name[0] == "#": continue
        fullname = os.path.join(init_path, name)
        basic_tools.doHOLA(fullname, config=config)

if __name__ == '__main__':
    run_tests('special/test_suite1')
    #run_tests('sbgn/test_suite1')
    #run_tests('metro/test_suite2')
