#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=80.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=80 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=80.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=85.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=85 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=85.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=90.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=90 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=90.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=95.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=95 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=95.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=100.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=100 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N5_k=100.txt

