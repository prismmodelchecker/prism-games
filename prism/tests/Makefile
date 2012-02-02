# PRISM executable (override at command-line to change)
PRISM_EXEC = prism

# Use TESTALL=true at command-line to not exit on first test fail
TESTALL =

MODEL_FILES = $(shell find . \( -name '*.pm' -o -name '*.sm' -o -name '*.nm' -o -name '*.smg' \) )
LOG_FILES = $(MODEL_FILES:%=%.log)

default: all

all: prism

prism: 
	time $(MAKE) prism_logs

prism_logs: $(LOG_FILES)

%.log: %
	@if [ -e $<.args ]; then \
	  MODEL_ARGS=`cat $<.args`; \
	fi; \
	if [ "$(TESTALL)" != "" ]; then \
	  TEST_ARG="-testall"; \
	else \
	  TEST_ARG="-test"; \
	fi; \
	NO_PROPS_EXIST=1; for PROP in $<*.props; do test -e $$PROP; NO_PROPS_EXIST=$$?; break; done; \
	if [ $$NO_PROPS_EXIST -eq 0 ]; then \
	  for PROP in $<*.props; do \
	    if [ -e $$PROP.args ]; then \
	      PROP_ARGS=`cat $$PROP.args`; \
	    fi; \
	    echo "$(PRISM_EXEC) $$TEST_ARG $$MODEL_ARGS $< $$PROP $$PROP_ARGS >> $@"; \
	    $(PRISM_EXEC) $$TEST_ARG $$MODEL_ARGS $< $$PROP $$PROP_ARGS >> $@ 2>&1; \
	    if [ $$? -ne 0 -a "$(TESTALL)" = "" ]; then \
	      cat $@; \
	      exit 1; \
	    fi; \
	    if [ "$(TESTALL)" != "" ]; then \
	     grep -E -i "error|fail" $@ || exit 0; \
	    fi; \
	  done; \
	else \
	  echo "$(PRISM_EXEC) $$TEST_ARG $$MODEL_ARGS $< > $@"; \
	  $(PRISM_EXEC) $$TEST_ARG $$MODEL_ARGS $< > $@ 2>&1; \
	  if [ $$? -ne 0 -a "$(TESTALL)" = "" ]; then \
	    cat $@; \
	    exit 1; \
	  fi; \
	  if [ "$(TESTALL)" != "" ]; then \
	     grep -E -i "error|fail" $@ || exit 0; \
	  fi; \
	fi \

clean:
	find . -name '*.log' -exec rm {} \;

celan: clean

#################################################
