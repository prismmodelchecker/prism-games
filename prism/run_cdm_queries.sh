#!/bin/sh

rm -f examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt

PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 > examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt
PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1,2>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 >> examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt
PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1,2,3>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 >> examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt
PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1,2,3,4>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 >> examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt
PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1,2,3,4,5>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 >> examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt
PRISM_JAVAMAXMEM=16g bin/prism examples/games/CDM/experiments/prob_choose_best/cdm7032.smg -ex -pctl '<<0,1,2,3,4,5,6>> Pmax=? [F all_committed&all_prefer_1]' -const Pexp=0.5,gamma=2.0,eta=2.0,lambda=2.0,Q1=1.0,Q2=0.5,Q3=0.25 >> examples/games/CDM/experiments/prob_choose_best/results/res_from0_g=2.0.txt


