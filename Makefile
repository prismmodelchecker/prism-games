############################################################
#  Small makefile for building PRISM source distributions  #
############################################################

default: none

none:
	@echo 'Did you want to build PRISM? Do "cd prism" and then "make"'

VERSION=# default value for VERSION is blank to force provision at command-line 

dist_src: dist_check_version
	mkdir dontcopy
	@if [ -e prism/examples ]; then \
	  echo "mv prism/examples dontcopy"; mv prism/examples dontcopy; \
	fi
	@if [ -e prism/tests ]; then \
	  echo "mv prism/tests dontcopy"; mv prism/tests dontcopy; \
	fi
	@if [ -e prism/examples-distr ]; then \
	  echo "mv prism/examples-distr prism/examples"; mv prism/examples-distr prism/examples; \
	else \
	  echo "mv prism-examples prism/examples"; mv prism-examples prism/examples; \
	fi
	mv cudd prism
	mv prism "prism-$(VERSION)-src"
	(cd "prism-$(VERSION)-src"; $(MAKE) dist_src )
	tar cfz "prism-$(VERSION)-src.tar.gz" "prism-$(VERSION)-src"

dist_check_version:
	@if [ "$(VERSION)" = "" ]; then echo "Usage: make dist_src VERSION=XXX"; exit 1; fi

#################################################
