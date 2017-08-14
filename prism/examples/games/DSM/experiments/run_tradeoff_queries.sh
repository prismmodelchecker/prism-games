#!/bin/sh




../../../../bin/prism value_tradeoff/dsm2304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res2.txt
../../../../bin/prism value_tradeoff/dsm2304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res2.txt




../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res3.txt
../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res3.txt
../../../../bin/prism value_tradeoff/dsm3304.smg -ex -pctl '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res3.txt




../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt
../../../../bin/prism value_tradeoff/dsm4304.smg -ex -pctl '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res4.txt




PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt
PRISM_JAVAMAXMEM=8g ../../../../bin/prism value_tradeoff/dsm5304.smg -ex -pctl '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]' >> value_tradeoff/results/res5.txt




PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt
PRISM_JAVAMAXMEM=16g ../../../../bin/prism value_tradeoff/dsm6304.smg -ex -pctl '<<p1,p2,p3,p4,p5,p6>> R{"value123456"}max=? [F time=max_time]'  >> value_tradeoff/results/res6.txt




PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1>> R{"value1"}max=? [F time=max_time]' > value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2>> R{"value12"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2,p3,p4,p5,p6>> R{"value123456"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt
PRISM_JAVAMAXMEM=32g ../../../../bin/prism value_tradeoff/dsm7304.smg -ex -pctl '<<p1,p2,p3,p4,p5,p6,p7>> R{"value1234567"}max=? [F time=max_time]' >> value_tradeoff/results/res7.txt


