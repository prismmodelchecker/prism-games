#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=5.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=5 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=5.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=10.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=10 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=10.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=15.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=15 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=15.txt

