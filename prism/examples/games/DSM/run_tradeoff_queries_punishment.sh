#!/bin/sh



../../../../bin/prism value_tradeoff/dsm2314p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres2.txt



../../../../bin/prism value_tradeoff/dsm3324p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres3.txt
../../../../bin/prism value_tradeoff/dsm3314p.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres3.txt



../../../../bin/prism value_tradeoff/dsm4334p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres4.txt
../../../../bin/prism value_tradeoff/dsm4324p.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres4.txt
../../../../bin/prism value_tradeoff/dsm4314p.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres4.txt



PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5344p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5334p.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5324p.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5314p.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres5.txt



PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6354p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6344p.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6334p.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6324p.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6314p.smg -ex -pctl '<<1,2,3,4,5>> R{"value12345"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres6.txt



PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7364p.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > examples/games/DSM/experiments/value_tradeoff/results/pres7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7354p.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7344p.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7334p.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7324p.smg -ex -pctl '<<1,2,3,4,5>> R{"value12345"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7314p.smg -ex -pctl '<<1,2,3,4,5,6>> R{"value123456"}max=? [F time=max_time]' >> examples/games/DSM/experiments/value_tradeoff/results/pres7.txt


