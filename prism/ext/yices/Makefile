################################################
#  NB: This Makefile is designed to be called  #
#      from the main PRISM Makefile. It won't  #
#      work on its own because it needs        #
#      various options to be passed in         #
################################################

default: all

all: checks yices

# MacOS libs
ifeq ($(OSTYPE),darwin)
    YICES_LIB = $(LIBPREFIX)yices.2$(LIBSUFFIX)
    YICES_JAVA_LIB = $(LIBPREFIX)yices2java$(LIBSUFFIX)
    GMP_LIB = $(LIBPREFIX)gmp.10$(LIBSUFFIX)
# Anything else (including unsupported OSs) - Linux libs
else
    YICES_LIB = $(LIBPREFIX)yices$(LIBSUFFIX).2.6
    YICES_JAVA_LIB = $(LIBPREFIX)yices2java$(LIBSUFFIX)
    GMP_LIB = $(LIBPREFIX)gmp$(LIBSUFFIX).10
endif
YICES_JAVA_JAR = yices.jar

# Try and prevent accidental makes (i.e. called manually, not from top-level Makefile)
checks:
	@if [ "$(LIBSUFFIX)" = "" ]; then \
	  (echo "Error: This Makefile is designed to be called from the main PRISM Makefile"; exit 1) \
	fi; 

yices:
	@(if [ "$(OSTYPE)" = "linux" ]; then \
	  echo "Copying Yices libraries for Linux"; \
	  cp "$(YICES_LIB)" "$(YICES_JAVA_LIB)" "$(YICES_JAVA_JAR)" "$(GMP_LIB)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)
	@(if [ "$(OSTYPE)" = "darwin" ]; then \
	  echo "Copying Yices libraries for MacOS"; \
	  cp "$(YICES_LIB)" "$(YICES_JAVA_LIB)" "$(YICES_JAVA_JAR)" "$(GMP_LIB)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)
	@(if [ "$(OSTYPE)" = "cygwin" ]; then \
	  echo "No support for Yices on Windows/Cywgin currently"; \
	  cp "$(YICES_JAVA_JAR)" ../../$(PRISM_LIB_DIR)/ ; \
	fi)

clean: checks
	rm -f "../../$(PRISM_LIB_DIR)/$(YICES_LIB)" "../../$(PRISM_LIB_DIR)/$(YICES_JAVA_LIB)" "../../$(PRISM_LIB_DIR)/$(YICES_JAVA_JAR)"

celan:	clean


#################################################
