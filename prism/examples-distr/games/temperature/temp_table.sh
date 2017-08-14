#!/bin/bash

export GREP_COLORS="ms=0;37:mc=01;32:sl=01;32:cx=01;31:fn=35:ln=:bn=32:se=36"

PROPS="phi_1 phi_2 phi_3 psi_1 psi_2 psi_3"
FILES="temperature.prism temperature.props temperature.sh"
echo -e "\e[1mTEMPERATURE CASE STUDY TIMING RESULTS\e[0m"
for f in $PROPS
do
	echo -e "\e[1m----\e[0m"
	echo -e "\e[1m$f:\e[0m"
	grep "Command line" temp_$f.log
	grep "Property constants" temp_$f.log
	grep --color=always -n "Synthesis took" temp_$f.log | tee /dev/tty | wc -l | \
	(
		read lines;
		if [ "$lines" -ne 3 ]
		then
			echo -e "\e[31mNOT SATISFIED\e[0m"
		fi
	)
	FILES="$FILES temp_$f.log temp_$f.strat"
done
echo -e "\e[1m----\e[0m"

# pack results in archive
tar -zcvf temp.tar.gz $FILES 2> /dev/null
echo -e "Results exported to \e[1mtemp.tar.gz\e[0m"
