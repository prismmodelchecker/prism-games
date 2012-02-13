#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=60.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=60 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=60.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=65.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=65 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=65.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=70.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=70 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=70.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=75.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=75 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=75.txt

