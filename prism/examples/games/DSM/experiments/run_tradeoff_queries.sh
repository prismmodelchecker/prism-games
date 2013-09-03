#!/bin/sh




../../../../bin/prism value_tradeoff/dsm2304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res2.txt
../../../../bin/prism value_tradeoff/dsm2304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res2.txt




../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res3.txt
../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res3.txt
../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res3.txt




../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt




PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<1,2,3,4,5>> R{"value12345"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt




PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1,2,3,4,5>> R{"value12345"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<1,2,3,4,5,6>> R{"value123456"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt




PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2,3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2,3,4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2,3,4,5>> R{"value12345"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2,3,4,5,6>> R{"value123456"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<1,2,3,4,5,6,7>> R{"value1234567"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt


