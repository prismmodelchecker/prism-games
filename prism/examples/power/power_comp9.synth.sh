LOG_ID=10

CITER=1000

# weight vector start index
DITEROFFSET1=1
DITEROFFSET2=1
DITEROFFSET3=1
DITEROFFSET4=1
DITEROFFSET5=1
DITEROFFSET6=1
DITEROFFSET7=1
DITEROFFSET8=1
DITEROFFSET9=1
DITEROFFSET10=1
DITEROFFSET11=1
DITEROFFSET12=1

# under which ID is strategy saved
STRAT_ID1=1
STRAT_ID2=2
STRAT_ID3=3
STRAT_ID4=4
STRAT_ID5=5
STRAT_ID6=6
STRAT_ID7=7
STRAT_ID8=8
STRAT_ID9=9
STRAT_ID10=10
STRAT_ID11=11
STRAT_ID12=12

# convergence epsilon
ACC1=0.01
ACC2=0.01
ACC3=0.001
ACC4=0.001
ACC5=0.1
ACC6=0.1
ACC7=0.1
ACC8=0.1
ACC9=0.1
ACC10=0.1
ACC11=0.1
ACC12=0.1

# if rounding used, the baseline accuracy
BACC1=100
BACC2=200
BACC3=1000
BACC4=2000
BACC5=100
BACC6=1000
BACC7=10000
BACC8=10000
BACC9=10000
BACC10=10000
BACC11=10000
BACC12=10000

# property index
PROP=5 # 5 for left, 9 for right

./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -logcpareto \
  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC1)) \
  -nocompatibility \
  -exportstrat $(($LOG_ID))_$(($STRAT_ID1)).strat \
  -multiditeroffset $DITEROFFSET1  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC1 \
  -prop $PROP \
  -gs 2>&1 1>pc9_$(($LOG_ID))_0.synth_$(($STRAT_ID1)).log &
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -logcpareto \
  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC2)) \
  -nocompatibility \
  -exportstrat $(($LOG_ID))_$(($STRAT_ID2)).strat \
  -multiditeroffset $DITEROFFSET2  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC2 \
  -prop $PROP \
  -gs 2>&1 1>pc9_$(($LOG_ID))_0.synth_$(($STRAT_ID2)).log &
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -logcpareto \
  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC3)) \
  -nocompatibility \
  -exportstrat $(($LOG_ID))_$(($STRAT_ID3)).strat \
  -multiditeroffset $DITEROFFSET3  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC3 \
  -prop $PROP \
  -gs 2>&1 1>pc9_$(($LOG_ID))_0.synth_$(($STRAT_ID3)).log &
./bin/prism ../prism-examples/power/power_comp9.smg ../prism-examples/power/power_comp9.props \
  -logcpareto \
  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC4)) \
  -nocompatibility \
  -exportstrat $(($LOG_ID))_$(($STRAT_ID4)).strat \
  -multiditeroffset $DITEROFFSET4  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC4 \
  -prop $PROP \
  -gs 2>&1 1>pc9_$(($LOG_ID))_0.synth_$(($STRAT_ID4)).log &


#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding \
#  -increasefactor 1.0 \
#  -baselineaccuracy $((BACC5)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID5)).strat \
#  -multiditeroffset $DITEROFFSET5 \
#  -multimaxditer 1 \
#  -multimaxciter $CITER \
#  -paretoepsilon $ACC5 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID5)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding \
#  -increasefactor 1.0 \
#  -baselineaccuracy $((BACC6)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID6)).strat \
#  -multiditeroffset $DITEROFFSET6 \
#  -multimaxditer 1 \
#  -multimaxciter $CITER \
#  -paretoepsilon $ACC6 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID6)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding \
#  -increasefactor 1.0 \
#  -baselineaccuracy $((BACC7)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID7)).strat \
#  -multiditeroffset $DITEROFFSET7 \
#  -multimaxditer 1 \
#  -multimaxciter $CITER \
#  -paretoepsilon $ACC7 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID7)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC8)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID8)).strat \
#  -multiditeroffset $DITEROFFSET8  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC8 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID8)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC9)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID9)).strat \
#  -multiditeroffset $DITEROFFSET9  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC9 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID9)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC10)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID10)).strat \
#  -multiditeroffset $DITEROFFSET10  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC10 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID10)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC11)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID11)).strat \
#  -multiditeroffset $DITEROFFSET11  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC11 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID11)).log &
#./bin/prism ../prism-examples/power/power_comp8.smg ../prism-examples/power/power_comp8.synth.props \
#  -logcpareto \
#  -multirounding  -increasefactor 1.0  -baselineaccuracy $((BACC12)) \
#  -nocompatibility \
#  -exportstrat $(($LOG_ID))_$(($STRAT_ID12)).strat \
#  -multiditeroffset $DITEROFFSET12  -multimaxditer 1  -multimaxciter $CITER  -paretoepsilon $ACC12 \
#  -prop $PROP \
#  -gs 2>&1 1>pc8_$(($LOG_ID))_0.synth_$(($STRAT_ID12)).log &
