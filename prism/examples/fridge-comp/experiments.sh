LOGFILE="exp"
LOGPOSTFIX=".log"
BASELINEACCURACY=200
INCREASEFACTOR=1.01

echo "Type in the log index, followed by [ENTER]:"
read LOGINDEX

echo "Type in the number of traders (1-4, or ALL for all traders), followed by [ENTER]:"
read TRADERS


if [[ $TRADERS = "1" || $TRADERS = "ALL" ]]; then

echo "Starting experiments with index "$LOGINDEX"."
echo ""
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 1 TRADER - COMPOSITIONAL                                  --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 1 trader - compositional
./bin/prism \
    -property 1 \
    -multirounding \
    -gs \
    -exportstrat rf1c.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-1.smg \
    ../prism-examples/fridge-comp/reach_fridge-1.props \
2>&1 | tee ${LOGFILE}1c${LOGINDEX}${LOGPOSTFIX}


echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 1 TRADER  - MONOLITHIC                                     --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 1 trader - monolithic
# note: no multiplier
./bin/prism \
    -property 2 \
    -multirounding \
    -gs \
    -exportstrat rf1m.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-1.smg \
    ../prism-examples/fridge-comp/reach_fridge-1.props \
2>&1 | tee ${LOGFILE}1m${LOGINDEX}${LOGPOSTFIX}

else

echo "Wrong number of traders entered."

fi


if [[ $TRADERS = "2" || $TRADERS = "ALL" ]]; then

echo ""
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 2 TRADERS - COMPOSITIONAL                                 --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 2 traders - compositional
./bin/prism \
    -property 1 \
    -multirounding \
    -gs \
    -exportstrat rf2c.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-2.smg \
    ../prism-examples/fridge-comp/reach_fridge-2.props \
2>&1 | tee ${LOGFILE}2c${LOGINDEX}${LOGPOSTFIX}


echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 2 TRADERS - MONOLITHIC                                    --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 2 traders - monolithic
# note: no multiplier
./bin/prism \
    -property 2 \
    -multirounding \
    -gs \
    -exportstrat rf2m.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-2.smg \
    ../prism-examples/fridge-comp/reach_fridge-2.props \
2>&1 | tee ${LOGFILE}2m${LOGINDEX}${LOGPOSTFIX}

else

echo "Wrong number of traders entered."

fi


if [[ $TRADERS = "3" || $TRADERS = "ALL" ]]; then

echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 3 TRADERS - COMPOSITIONAL                                 --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 3 traders - compositional
./bin/prism \
    -property 1 \
    -multirounding \
    -gs \
    -exportstrat rf3c.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-3.smg \
    ../prism-examples/fridge-comp/reach_fridge-3.props \
2>&1 | tee ${LOGFILE}3c${LOGINDEX}${LOGPOSTFIX}


echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 3 TRADERS - MONOLITHIC                                    --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 3 traders - monolithic
./bin/prism \
    -property 2 \
    -multirounding \
    -gs \
    -exportstrat rf3m.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-3.smg \
    ../prism-examples/fridge-comp/reach_fridge-3.props \
2>&1 | tee ${LOGFILE}3m${LOGINDEX}${LOGPOSTFIX}

else

echo "Wrong number of traders entered."

fi


if [[ $TRADERS = "4" || $TRADERS = "ALL" ]]; then

echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 4 TRADERS - COMPOSITIONAL                                 --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 4 traders - compositional
./bin/prism \
    -property 1 \
    -multirounding \
    -gs \
    -exportstrat rf4c.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-4.smg \
    ../prism-examples/fridge-comp/reach_fridge-4.props \
2>&1 | tee ${LOGFILE}4c${LOGINDEX}${LOGPOSTFIX}


echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"
echo "--                                                           --"
echo "-- 4 TRADERS - MONOLITHIC                                    --"
echo "--                                                           --"
echo "---------------------------------------------------------------"
echo "---------------------------------------------------------------"

# 4 traders - monolithic
./bin/prism \
    -property 2 \
    -multirounding \
    -gs \
    -exportstrat rf4m.strat \
    -baselineaccuracy $BASELINEACCURACY \
    -increasefactor $INCREASEFACTOR \
    ../prism-examples/fridge-comp/reach_fridge-4.smg \
    ../prism-examples/fridge-comp/reach_fridge-4.props \
2>&1 | tee ${LOGFILE}4m${LOGINDEX}${LOGPOSTFIX}

else

echo "Wrong number of traders entered."

fi
