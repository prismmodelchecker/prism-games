#!/bin/bash

export GREP_COLORS="ms=0;37:mc=01;37:sl=01;32:cx=0;37:fn=35:ln=:bn=32:se=36"

echo -e "\e[1mPOWER CASE STUDY TIMING RESULTS\e[0m"
PROPS="1 2 3 4 5 6 7 8 9 10 11 12"
FILES="power.prism power.props power.sh"
for f in $PROPS
do
	echo -e "\e[1m----\e[0m"
	echo -e "\e[1mpower_14_$f:\e[0m"
	grep "Command line" power_14_$f.log | \
	  grep -oh -E "N=[0-9]*|del_max=[0-9]*|baselineaccuracy [0-9]*|paretoepsilon [0-9]*.[0-9]*"
	grep "Property constants" power_14_$f.log
	grep --color "Computing reachable states..." power_14_$f.log
	grep -B 4 --color=always -n "Synthesis took" power_14_$f.log | \
          grep -E "C-ITER|Synthesis took" | \
	  tee /dev/tty | wc -l | \
	(
		read lines;
		if [ "$lines" -ne 4 ]
		then
			echo -e "\e[31mNOT SATISFIED\e[0m"
		fi
	)
	FILES="$FILES power_14_$f.log power_14_$f.strat"
done
echo -e "\e[1m----\e[0m"

# pack results in archive
tar -zcf power.tar.gz $FILES 2> /dev/null
echo -e "Results exported to \e[1mpower.tar.gz\e[0m"
