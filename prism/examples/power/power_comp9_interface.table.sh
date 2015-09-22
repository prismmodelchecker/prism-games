LOG_ID=10

CITER=200
DITER=100

# convergence epsilon
ACC1=0.01
ACC2=0.01
ACC3=0.01
ACC4=0.01
ACC5=0.01
ACC6=0.01

# if rounding used, the baseline accuracy
BACC1=100
BACC2=100
BACC3=100
BACC4=100
BACC5=100
BACC6=100

# property index
PROP_INTERFACE=9 # 3 for both, 9 for left only

# Experiment 7
EXP_ID=7
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.85,l_buses_r=0.85,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 8
EXP_ID=8
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.6,l_i1_r=0.6 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 9
EXP_ID=9
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.85,l_buses_r=0.85,l_i1_l=0.7,l_i1_r=0.7 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 10
EXP_ID=10
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.7,l_i1_r=0.7 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 11
EXP_ID=11
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.85,l_buses_r=0.85,l_i1_l=0.8,l_i1_r=0.8 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &

# Experiment 12
EXP_ID=12
./bin/prism ../prism-examples/power/power_comp9_interface.smg ../prism-examples/power/power_comp9.props \
  -prop $PROP_INTERFACE \
  -const max_gc=0,u_fail_l=0.02,u_fail_r=0.02,l_buses_l=0.9,l_buses_r=0.9,l_i1_l=0.8,l_i1_r=0.8 \
  -logcpareto -nocompatibility -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -multimaxciter $CITER -paretoepsilon $ACC1 -multimaxditer $DITER \
  -exportstrat pc9_$(($LOG_ID))_$(($EXP_ID)).strat -gs 2>&1 1>pc9_$(($LOG_ID))_$(($EXP_ID)).log &
