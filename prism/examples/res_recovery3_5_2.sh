#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=20.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=20 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=20.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=25.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=25 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=25.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=30.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=30 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=30.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=35.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=35 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=35.txt

