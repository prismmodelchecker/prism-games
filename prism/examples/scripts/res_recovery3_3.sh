#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=10.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=10 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=10.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=20.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=20 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=20.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=30.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=30 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=30.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=40.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=40 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=40.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=50.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=50 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=50.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=60.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=60 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=60.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=70.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=70 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=70.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=80.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=80 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=80.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=90.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=90 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=90.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=100.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm3032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N3.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=100 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N3_k=100.txt

