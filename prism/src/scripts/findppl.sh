#!/bin/sh

# Detection of PPL:
# * By default, try to use the ppl-config command (if it's in the path)
# * If that fails, try looking for an rpm

if command ppl-config &> /dev/null; then
  echo `ppl-config -l`
else
  if command rpm &> /dev/null; then
	echo `rpm -qi --filesbypkg ppl-java | grep "libppl_java" | awk -F'/' '{OFS="/";$1="";$(NF-1)="";$NF="";print}' | sed 's/.$//'`
  fi
fi
