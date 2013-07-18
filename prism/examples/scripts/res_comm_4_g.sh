#!/bin/sh

rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=0e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=0,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=0e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=1e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=1,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=1e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=2e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=2,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=2e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=3e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=3,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=3e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=4e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=4,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=4e=1.txt
rm -f examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=5e=1.txt
PRISM_JAVAMAXMEM=16g bin/prism -ex examples/games/CDM/experiments/no_commit/cdm4032.smg examples/games/CDM/experiments/no_commit/props_comm_4.pctl -const Pexp=0.5,Q1=1.0,Q2=0.5,Q3=0.25,gamma=5,eta=1,lambda=1 > examples/games/CDM/experiments/no_commit/results/res_comm_4_l=1g=5e=1.txt

