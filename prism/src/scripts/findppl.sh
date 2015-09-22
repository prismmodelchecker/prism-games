#!/bin/sh

if type rpm &> /dev/null; then
  echo `rpm -qi --filesbypkg ppl-java | grep "libppl_java" | awk -F'/' '{OFS="/";$1="";$(NF-1)="";$NF="";print}' | sed 's/.$//'`
else
  echo `ppl-config -p`/ppl/lib
fi
