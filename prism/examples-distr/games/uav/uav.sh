#!/bin/bash

TIMEOUT=240 # in minutes

if hash timeout 2>/dev/null; then
        TO=timeout
else
        TO=gtimeout
fi

$TO $((TIMEOUT))m \
../../bin/prism uav{.prism,.props} -prop 2 -logcpareto -const accu_load1=0.9,accu_load2=0.8,fd=0.7,COUNTER=2,del=0.5 -multirounding -increasefactor 1.0 -baselineaccuracy 1000 -multimaxciter 100 -paretoepsilon 0.01 -pareto -gs 2>&1 | tee uav.log

