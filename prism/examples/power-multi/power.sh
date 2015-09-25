#!/bin/bash

LOG_ID=14
TIMEOUT=240 # in minutes

if hash timeout 2>/dev/null; then
        TO=timeout
else
        TO=gtimeout
fi

CITER=1000
MAXM=100000

# convergence epsilon
ACC1=0.001
ACC2=0.001
ACC3=0.001
ACC4=0.001
ACC5=0.001
ACC6=0.001
ACC7=0.01
ACC8=0.01
ACC9=0.01
ACC10=0.01
ACC11=0.01
ACC12=0.01

# if rounding used, the baseline accuracy
BACC1=1000
BACC2=1000
BACC3=1000
BACC4=1000
BACC5=1000
BACC6=1000
BACC7=100
BACC8=100
BACC9=100
BACC10=100
BACC11=100
BACC12=100

# property index
PROP_NO_INTERFACE=1
PROP_INTERFACE=2

{

# Experiment 1
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=0,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC1 \
  -exportstrat power_$(($LOG_ID))_1.strat -gs 2>&1 1> power_$(($LOG_ID))_1.log ;

# Experiment 2
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=1,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC2)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC2 \
  -exportstrat power_$(($LOG_ID))_2.strat -gs 2>&1 1> power_$(($LOG_ID))_2.log ;

# Experiment 3
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=1,del_max=1,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC3)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC3 \
  -exportstrat power_$(($LOG_ID))_3.strat -gs 2>&1 1> power_$(($LOG_ID))_3.log ;

# Experiment 4
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=2,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC4)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC4 \
  -exportstrat power_$(($LOG_ID))_4.strat -gs 2>&1 1> power_$(($LOG_ID))_4.log ;

# Experiment 5
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=2,del_max=1,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC5)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC5 \
  -exportstrat power_$(($LOG_ID))_5.strat -gs 2>&1 1> power_$(($LOG_ID))_5.log ;

# Experiment 6
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,N=2,del_max=2,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC6)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC6 \
  -exportstrat power_$(($LOG_ID))_6.strat -gs 2>&1 1> power_$(($LOG_ID))_6.log ;

# Experiment 7
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=0,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC7)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC7 \
  -exportstrat power_$(($LOG_ID))_7.strat -gs 2>&1 1> power_$(($LOG_ID))_7.log ;

# Experiment 8
timeout $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=1,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.95,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC8)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC8 \
  -exportstrat power_$(($LOG_ID))_8.strat -gs 2>&1 1> power_$(($LOG_ID))_8.log ;

# Experiment 9
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=1,del_max=1,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC9)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC9 \
  -exportstrat power_$(($LOG_ID))_9.strat -gs 2>&1 1> power_$(($LOG_ID))_9.log ;

# Experiment 10
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=2,del_max=0,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC10)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC10 \
  -exportstrat power_$(($LOG_ID))_10.strat -gs 2>&1 1> power_$(($LOG_ID))_10.log ; 

# Experiment 11
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=2,del_max=1,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC11)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC11 \
  -exportstrat power_$(($LOG_ID))_11.strat -gs 2>&1 1> power_$(($LOG_ID))_11.log ;

# Experiment 12
$TO $((TIMEOUT))m \
../../bin/prism power{.prism,.props} \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,N=2,del_max=2,u_fail_l=0.01,u_fail_r=0.01,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC12)) \
  -multimaxciter $CITER -multimaxm $MAXM -paretoepsilon $ACC12 \
  -exportstrat power_$(($LOG_ID))_12.strat -gs 2>&1 1> power_$(($LOG_ID))_12.log 

} &
