#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=0g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=0 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=0g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=1g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=1g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=2g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=2 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=2g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=3g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=3 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=3g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=4g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=4 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=4g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_total_4_l=5g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_total_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=5 > examples/games/CDM/experiments/no_commit/results/res_total_4_l=5g=1e=1.txt


