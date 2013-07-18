#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=40.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=40 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=40.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=45.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=45 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=45.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=50.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=50 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=50.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=55.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=55 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=55.txt

