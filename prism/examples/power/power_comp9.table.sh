LOG_ID=10

CITER=1000

# convergence epsilon
ACC1=0.001
ACC2=0.001
ACC3=0.001
ACC4=0.001
ACC5=0.01
ACC6=0.01
ACC7=0.01
ACC8=0.01

# if rounding used, the baseline accuracy
BACC1=1000
BACC2=1000
BACC3=1000
BACC4=1000
BACC5=100
BACC6=100
BACC7=100
BACC8=100

# property index
PROP_NO_INTERFACE=7 # 1 for both, 7 for left only
PROP_INTERFACE=8 # 2 for both, 8 for left only

# Experiment 1
EXP_ID=1
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,max_gc=0,max_delay=0,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.85,l_buses_r=0.85,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 2
EXP_ID=2
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,max_gc=0,max_delay=0,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC2)) \
  -multimaxciter $CITER -paretoepsilon $ACC2 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 3
EXP_ID=3
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,max_gc=3,max_delay=1,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.85,l_buses_r=0.85,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC3)) \
  -multimaxciter $CITER -paretoepsilon $ACC3 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 4
EXP_ID=4
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_NO_INTERFACE \
  -const I1_health=0.0,max_gc=3,max_delay=1,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC4)) \
  -multimaxciter $CITER -paretoepsilon $ACC4 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 5
EXP_ID=5
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,max_gc=0,max_delay=0,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC5)) \
  -multimaxciter $CITER -paretoepsilon $ACC5 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 6
EXP_ID=6
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,max_gc=0,max_delay=0,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.95,l_buses_r=0.95,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC6)) \
  -multimaxciter $CITER -paretoepsilon $ACC6 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 7
EXP_ID=7
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,max_gc=2,max_delay=1,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC7)) \
  -multimaxciter $CITER -paretoepsilon $ACC7 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 8
EXP_ID=8
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const I1_health=0.6,max_gc=2,max_delay=1,u_fail_l=0.01,u_fail_r=0.02,l_buses_l=0.95,l_buses_r=0.95,l_i1_l=0.0,l_i1_r=0.0 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC8)) \
  -multimaxciter $CITER -paretoepsilon $ACC8 \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

