#!/bin/bash

# probability at least half sent
../../bin/prism jamming.prism jamming.props -javamaxmem 8g -prop 1 -const slots=5:5:25
# expected number sent
../../bin/prism jamming.prism jamming.props -javamaxmem 8g -prop 2 -const slots=5:5:25
