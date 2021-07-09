################################################
#  NB: This Makefile is designed to be called  #
#      from the main PRISM Makefile. It won't  #
#      work on its own because it needs        #
#      various options to be passed in         #
################################################

default: all

all: checks z3

# Cygwin libs
ifeq ($(OSTYPE),cygwin)
    Z3_LIB = $(LIBPREFIX)libz3$(LIBSUFFIX)
    Z3_JAVA_LIB = $(LIBPREFIX)libz3java$(LIBSUFFIX)
# Anything else (including unsupported OSs) - Mac/Linux libs
else
    Z3_LIB = $(LIBPREFIX)z3$(LIBSUFFIX)
    Z3_JAVA_LIB = $(LIBPREFIX)z3java$(LIBSUFFIX)
endif
Z3_JAVA_JAR = com.microsoft.z3.jar

# Try and prevent accidental makes (i.e. called manually, not from top-level Makefile)
checks:
	@if [ "$(LIBSUFFIX)" = "" ]; then \
	  (echo "Error: This Makefile is designed to be called from the main PRISM Makefile"; exit 1) \
	fi; 

z3:
	@(if [ "$(OSTYPE)" = "linux" ]; then \
	  echo "Copying Z3 libraries for Linux"; \
	  cp "$(Z3_LIB)" "$(Z3_JAVA_LIB)" "$(Z3_JAVA_JAR)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)
	@(if [ "$(OSTYPE)" = "darwin" ]; then \
	  echo "Copying Z3 libraries for MacOS"; \
	  cp "$(Z3_LIB)" "$(Z3_JAVA_LIB)" "$(Z3_JAVA_JAR)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)
	@(if [ "$(OSTYPE)" = "cygwin" ]; then \
	  echo "Copying Z3 libraries for Cywgin"; \
	  cp "$(Z3_LIB)" "$(Z3_JAVA_LIB)" "$(Z3_JAVA_JAR)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)

clean: checks
	rm -f "../../$(PRISM_LIB_DIR)/$(Z3_LIB)" "../../$(PRISM_LIB_DIR)/$(Z3_JAVA_LIB)" "../../$(PRISM_LIB_DIR)/$(Z3_JAVA_JAR)"

celan:	clean


#################################################
