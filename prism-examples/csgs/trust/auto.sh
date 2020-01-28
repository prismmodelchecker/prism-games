#!/bin/bash

# expected number of unpaid

../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 1 -const td=0,K=1:1:9
../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 1 -const td=1,K=1:1:9
../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 1 -const td=2,K=1:1:9

# radio of unpaid

../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 2 -const td=0,K=1:1:9
../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 2 -const td=1,K=1:1:9
../../bin/prism trust_csg.prism trust_csg.props -javamaxmem 8g -prop 2 -const td=2,K=1:1:9
