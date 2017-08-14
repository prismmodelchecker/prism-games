#!/bin/bash

# assume-guarantee synthesis
../../bin/prism comp_example{.prism,.props} -prop 1 \
     -multimaxciter 500 -multimaxditer 500 -paretoepsilon 0.001 -logcpareto -logdpareto -gs \
     -exportstrat comp_example.strat \
    2>&1 | tee comp_example.log

# monolithic synthesis
../../bin/prism comp_example{.prism,.props} -prop 4 \
     -multimaxciter 500 -multimaxditer 500 -paretoepsilon 0.001 -logcpareto -gs \
     -exportstrat comp_example_mono.strat \
    2>&1 | tee -a comp_example.log
