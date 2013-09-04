#!/bin/sh

prism dsm2314p.smg '<<p1>> R{"value1"}max=? [F time=max_time]' 
prism dsm2304p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]'

prism dsm3324p.smg '<<p1>> R{"value1"}max=? [F time=max_time]' 
prism dsm3314p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]' 
prism dsm3304p.smg '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' 

prism dsm4334p.smg '<<p1>> R{"value1"}max=? [F time=max_time]' 
prism dsm4324p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]' 
prism dsm4314p.smg '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' 
prism dsm4304p.smg '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]' 

prism dsm5344p.smg '<<p1>> R{"value1"}max=? [F time=max_time]' 
prism dsm5334p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]' 
prism dsm5324p.smg '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' 
prism dsm5314p.smg '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]'
prism dsm5304p.smg '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]'

prism dsm6354p.smg '<<p1>> R{"value1"}max=? [F time=max_time]'
prism dsm6344p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]'
prism dsm6334p.smg '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]'
prism dsm6324p.smg '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]'
prism dsm6314p.smg '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]'
prism dsm6304p.smg '<<p1,p2,p3,p4,p5,p6>> R{"value123456"}max=? [F time=max_time]'

prism dsm7364p.smg '<<p1>> R{"value1"}max=? [F time=max_time]'
prism dsm7354p.smg '<<p1,p2>> R{"value12"}max=? [F time=max_time]'
prism dsm7344p.smg '<<p1,p2,p3>> R{"value123"}max=? [F time=max_time]' 
prism dsm7334p.smg '<<p1,p2,p3,p4>> R{"value1234"}max=? [F time=max_time]' 
prism dsm7324p.smg '<<p1,p2,p3,p4,p5>> R{"value12345"}max=? [F time=max_time]' 
prism dsm7314p.smg '<<p1,p2,p3,p4,p5,p6>> R{"value123456"}max=? [F time=max_time]' 
prism dsm7304p.smg '<<p1,p2,p3,p4,p5,p6,p7>> R{"value1234567"}max=? [F time=max_time]' 

