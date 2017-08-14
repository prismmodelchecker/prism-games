#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=5.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=5 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=5.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=10.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=10 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=10.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=15.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=15 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=15.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=20.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=20 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=20.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=25.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=25 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=25.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=30.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=30 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=30.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=35.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=35 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=35.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=40.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=40 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=40.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=45.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=45 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=45.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=50.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=50 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=50.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=55.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=55 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=55.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=60.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=60 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=60.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=65.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=65 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=65.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=70.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=70 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=70.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=75.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=75 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=75.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=80.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=80 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=80.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=85.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=85 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=85.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=90.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=90 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=90.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=95.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=95 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=95.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=100.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_3_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=100 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_3_N4_k=100.txt

