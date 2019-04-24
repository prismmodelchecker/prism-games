#!/bin/bash
TIMEOUT=240

if hash timeout 2>/dev/null; then
	TO=timeout
else
	TO=gtimeout
fi

{
$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 7 \
	-multimaxciter 100 -multimaxditer 20 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 100 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_phi_1.strat \
	-paretoepsilon 0.05 \
	2>&1 1> temp_phi_1.log ; 

$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 7 \
	-multimaxciter 100 -multimaxditer 20 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 500 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_phi_2.strat \
	-paretoepsilon 0.02 \
	2>&1 1> temp_phi_2.log ; 

$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 7 \
	-multimaxciter 100 -multimaxditer 20 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 1000 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_phi_3.strat \
	-paretoepsilon 0.01 \
	2>&1 1> temp_phi_3.log ; 

$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 8 \
	-multimaxciter 250 -multimaxditer 10 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 100 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_psi_1.strat \
	-paretoepsilon 0.05 \
	2>&1 1> temp_psi_1.log ; 

$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 8 \
	-multimaxciter 250 -multimaxditer 10 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 500 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_psi_2.strat \
	-paretoepsilon 0.02 \
	2>&1 1> temp_psi_2.log ; 

$TO $((TIMEOUT))m \
../../bin/prism temperature{.prism,.props} \
	-prop 8 \
	-multimaxciter 250 -multimaxditer 10 -gs \
	-multiminm 2 -multimaxm 10000 \
	-baselineaccuracy 1000 -increasefactor 1.0 -multirounding \
	-logcpareto -logdpareto \
	-exportstrat temp_psi_3.strat \
	-paretoepsilon 0.01 \
	2>&1 1> temp_psi_3.log ; 

} &
