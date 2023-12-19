################################################
#  NB: This Makefile is designed to be called  #
#      from the main PRISM Makefile. It won't  #
#      work on its own because it needs        #
#      various options to be passed in         #
################################################

default: all

all: checks z3

# Z3 library filenames
Z3_LIB = z3_lib_missing
Z3_JAVA_LIB = z3_java_lib_missing
ifeq ($(OSTYPE),linux)
    Z3_LIB = $(LIBPREFIX)z3$(LIBSUFFIX)
    Z3_JAVA_LIB = $(LIBPREFIX)z3java$(LIBSUFFIX)
endif
ifeq ($(OSTYPE),darwin)
    Z3_LIB = $(LIBPREFIX)z3$(LIBSUFFIX)
    Z3_JAVA_LIB = $(LIBPREFIX)z3java$(LIBSUFFIX)
endif
ifeq ($(OSTYPE),cygwin)
    Z3_LIB = $(LIBPREFIX)libz3$(LIBSUFFIX)
    Z3_JAVA_LIB = $(LIBPREFIX)libz3java$(LIBSUFFIX)
endif
Z3_JAVA_JAR = com.microsoft.z3.jar

# Try and prevent accidental makes (i.e. called manually, not from top-level Makefile)
checks:
	@if [ "$(LIBSUFFIX)" = "" ]; then \
	  (echo "Error: This Makefile is designed to be called from the main PRISM Makefile"; exit 1) \
	fi; 

# Copy files if they exist and are newer
z3: ../../$(PRISM_LIB_DIR)/$(Z3_LIB) ../../$(PRISM_LIB_DIR)/$(Z3_JAVA_LIB) ../../$(PRISM_LIB_DIR)/$(Z3_JAVA_JAR)

../../$(PRISM_LIB_DIR)/$(Z3_LIB): $(OSTYPE)/$(ARCH)/$(Z3_LIB)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

../../$(PRISM_LIB_DIR)/$(Z3_JAVA_LIB): $(OSTYPE)/$(ARCH)/$(Z3_JAVA_LIB)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

../../$(PRISM_LIB_DIR)/$(Z3_JAVA_JAR): $(Z3_JAVA_JAR)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

$(OSTYPE)/$(ARCH)/z3_lib_missing: missing

$(OSTYPE)/$(ARCH)/z3_java_lib_missing: missing

missing:
	@echo "Missing pre-built native libraries for $(OSTYPE) $(ARCH)"

clean: checks
	rm -f "../../$(PRISM_LIB_DIR)/$(Z3_LIB)"
	rm -f "../../$(PRISM_LIB_DIR)/$(Z3_JAVA_LIB)"
	rm -f "../../$(PRISM_LIB_DIR)/$(Z3_JAVA_JAR)"

celan:	clean


#################################################
