#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=1 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=2.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=2 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=2.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=3.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=3 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=3.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=4.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=4 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=4.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=5.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=5 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=5.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=6.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=6 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=6.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=7.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=7 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=7.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=8.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=8 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=8.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=9.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=9 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=9.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=10.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=10 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=10.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=11.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=11 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=11.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=12.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=12 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=12.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=13.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=13 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=13.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=14.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=14 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=14.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=15.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=15 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=15.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=16.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=16 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=16.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=17.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=17 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=17.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=18.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=18 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=18.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=19.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=19 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=19.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=20.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=20 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=20.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=21.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=21 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=21.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=22.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=22 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=22.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=23.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=23 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=23.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=24.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=24 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=24.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=25.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=25 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=25.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=26.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=26 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=26.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=27.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=27 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=27.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=28.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=28 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=28.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=29.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=29 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=29.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=30.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=30 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=30.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=31.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=31 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=31.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=32.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=32 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=32.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=33.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=33 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=33.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=34.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=34 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=34.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=35.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=35 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=35.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=36.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=36 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=36.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=37.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=37 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=37.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=38.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=38 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=38.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=39.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=39 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=39.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=40.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=40 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=40.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=41.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=41 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=41.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=42.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=42 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=42.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=43.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=43 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=43.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=44.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=44 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=44.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=45.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=45 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=45.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=46.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=46 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=46.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=47.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=47 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=47.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=48.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=48 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=48.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=49.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=49 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=49.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=50.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_recovery_from_2_N4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0,k=50 > examples/games/CDM/experiments/no_commit/results/res_recovery_from_2_N4_k=50.txt

