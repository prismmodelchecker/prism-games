#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=5.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=5 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=5.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=10.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=10 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=10.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=15.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=15 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=15.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=20.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=20 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=20.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=25.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=25 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=25.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=30.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=30 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=30.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=35.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=35 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=35.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=40.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=40 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=40.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=45.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=45 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=45.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=50.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=50 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=50.txt

