#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=0.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=0.0,eta=1.0,lambda=1.0,k=0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=0.0e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=0.5e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=0.5,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=0.5e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=1.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1.0,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=1.0e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=1.5e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1.5,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=1.5e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=2.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.0,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=2.0e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=2.5e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2.5,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=2.5e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=3.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=3.0,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=3.0e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=3.5e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=3.5,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=3.5e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=4.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=4.0,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=4.0e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=4.5e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=4.5,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=4.5e=1.0.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=5.0e=1.0.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm5032.smg examples/games/CDM/experiments/no_commit/props_runtime_5.pctl -const k=0,Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=5.0,eta=1.0,lambda=1.0 > examples/games/CDM/experiments/no_commit/results/res_runtime_5_l=1.0g=5.0e=1.0.txt

