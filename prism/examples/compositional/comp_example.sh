#!/bin/bash

# assume-guarantee synthesis
../../prism/bin/prism comp_example{.smg,.props} -prop 1 \
     -multimaxciter 500 -multimaxditer 500 -paretoepsilon 0.001 -logcpareto -gs \
     -exportstrat comp_example.strat \
    2>&1 | tee comp_example.log

# monolithic synthesis
../../prism/bin/prism comp_example{.smg,.props} -prop 4 \
     -multimaxciter 500 -multimaxditer 500 -paretoepsilon 0.001 -logcpareto -gs \
     -exportstrat comp_example.strat \
    2>&1 | tee -a comp_example.log
