################################################
#  NB: This Makefile is designed to be called  #
#      from the main PRISM Makefile. It won't  #
#      work on its own because it needs        #
#      various options to be passed in         #
################################################

default: all

all: checks yices

# Yices library filenames
YICES_LIB = yices_lib_missing
YICES_JAVA_LIB = yices_java_lib_missing
ifeq ($(OSTYPE),linux)
    YICES_LIB = $(LIBPREFIX)yices$(LIBSUFFIX).2.6
    YICES_JAVA_LIB = $(LIBPREFIX)yices2java$(LIBSUFFIX)
endif
ifeq ($(OSTYPE),darwin)
    YICES_LIB = $(LIBPREFIX)yices.2$(LIBSUFFIX)
    YICES_JAVA_LIB = $(LIBPREFIX)yices2java$(LIBSUFFIX)
endif
ifeq ($(OSTYPE),cygwin)
    YICES_LIB = $(LIBPREFIX)libyices$(LIBSUFFIX)
    YICES_JAVA_LIB = $(LIBPREFIX)yices2java$(LIBSUFFIX)
endif
YICES_JAVA_JAR = yices.jar

# Try and prevent accidental makes (i.e. called manually, not from top-level Makefile)
checks:
	@if [ "$(LIBSUFFIX)" = "" ]; then \
	  (echo "Error: This Makefile is designed to be called from the main PRISM Makefile"; exit 1) \
	fi; 

# Copy files if they exist and are newer
yices: ../../$(PRISM_LIB_DIR)/$(YICES_LIB) ../../$(PRISM_LIB_DIR)/$(YICES_JAVA_LIB) ../../$(PRISM_LIB_DIR)/$(YICES_JAVA_JAR)

../../$(PRISM_LIB_DIR)/$(YICES_LIB): $(OSTYPE)/$(ARCH)/$(YICES_LIB)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

../../$(PRISM_LIB_DIR)/$(YICES_JAVA_LIB): $(OSTYPE)/$(ARCH)/$(YICES_JAVA_LIB)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

../../$(PRISM_LIB_DIR)/$(YICES_JAVA_JAR): $(YICES_JAVA_JAR)
	@if [ -f "$<" ]; then \
	  echo cp $< $@; \
	  cp $< $@; \
	fi

$(OSTYPE)/$(ARCH)/yices_lib_missing: missing

$(OSTYPE)/$(ARCH)/yices_java_lib_missing: missing

missing:
	@echo "Missing pre-built native libraries for $(OSTYPE) $(ARCH)"

clean: checks
	rm -f "../../$(PRISM_LIB_DIR)/$(YICES_LIB)"
	rm -f "../../$(PRISM_LIB_DIR)/$(YICES_JAVA_LIB)"
	rm -f "../../$(PRISM_LIB_DIR)/$(YICES_JAVA_JAR)"

celan:	clean


#################################################
