#./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -multiditeroffset 0 -multimaxditer 1 -multimaxciter 200 -paretoepsilon 0.001 -prop 10 2>&1 | tee pc5_1.log

CITER=21
DITER=100
LOG_ID=28
ACC=0.01
PROP=13 # 7 for HVAC left, 13 for HVAC right

#./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -exportstrat 0.strat -multiditeroffset $(($DITER*0)) -multimaxditer $DITER -multimaxciter $CITER -paretoepsilon $ACC -prop $PROP 2>&1 1>pc5_$(($LOG_ID))_0.log &


./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -multiditeroffset $(($DITER*0)) -multimaxditer $DITER -multimaxciter $CITER -paretoepsilon $ACC -prop $PROP 2>&1 1>pc5_$(($LOG_ID))_0.log &
./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -multiditeroffset $(($DITER*1)) -multimaxditer $DITER -multimaxciter $CITER -paretoepsilon $ACC -prop $PROP 2>&1 1>pc5_$(($LOG_ID))_1.log &
./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -multiditeroffset $(($DITER*2)) -multimaxditer $DITER -multimaxciter $CITER -paretoepsilon $ACC -prop $PROP 2>&1 1>pc5_$(($LOG_ID))_2.log &
./bin/prism ../prism-examples/power/power_comp5.smg ../prism-examples/power/power_comp5.props -multiditeroffset $(($DITER*3)) -multimaxditer $DITER -multimaxciter $CITER -paretoepsilon $ACC -prop $PROP 2>&1 1>pc5_$(($LOG_ID))_3.log &

