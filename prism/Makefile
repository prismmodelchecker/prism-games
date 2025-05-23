##############################################
#  NB: This is the main Makefile for PRISM.  #
#      It calls all the other Makefiles in   #
#      subdirectories, passing in all the    #
#      options configured here.              #
##############################################

###############
# Directories #
###############

# Note that these are all relative to the PRISM directory
# to make the distribution more 'portable'.
# If this is a problem, the best solution is to create symlinks.

export PRISM_SRC_DIR     = src
export PRISM_CLASSES_DIR = classes
export PRISM_OBJ_DIR     = obj
export PRISM_LIB_DIR     = lib
export PRISM_INCLUDE_DIR = include
export PRISM_IMAGES_DIR  = images
export PRISM_DTDS_DIR    = dtds
export PRISM_TESTS_DIR   = unit-tests

# Location of CUDD (used to be variable; now mainly fixed with the git repo layout)

CUDD_DIR = ../cudd

####################
# Operating system #
####################

# OSTYPE needs to be one of: linux, solaris, cygwin, darwin
# This makefile will try to detect which one of these is appropriate.
# If this detection does not work, or you wish to override it,
# either uncomment one of the lines directly below
# or pass a value to make directly, e.g.: make OSTYPE=linux

#OSTYPE = linux
#OSTYPE = solaris
#OSTYPE = cygwin
#OSTYPE = darwin

ifdef OSTYPE
	# Look for common variants, e.g. gnu-linux -> linux
	ifneq (,$(findstring linux, $(OSTYPE)))
	  OSTYPE = linux
	endif
	ifneq (,$(findstring solaris, $(OSTYPE)))
	  OSTYPE = solaris
	endif
	ifneq (,$(findstring cygwin, $(OSTYPE)))
	  OSTYPE = cygwin
	endif
	# For Cygwin , OSTYPE is sometimes set to "posix"
	ifneq (,$(findstring posix, $(OSTYPE)))
	  OSTYPE = cygwin
	endif
	ifneq (,$(findstring darwin, $(OSTYPE)))
	  OSTYPE = darwin
	endif
else
	# If OSTYPE is not defined/available, try uname
	ifneq (,$(findstring Linux, $(shell uname -s)))
		OSTYPE = linux
	endif
	ifneq (,$(findstring SunOS, $(shell uname -s)))
		OSTYPE = solaris
	endif
	ifneq (,$(findstring CYGWIN, $(shell uname -s)))
		OSTYPE = cygwin
	endif
	ifneq (,$(findstring Darwin, $(shell uname -s)))
		OSTYPE = darwin
	endif
endif

################
# Architecture #
################

# For Linux/Mac, we use uname to check the architecture
ifeq ($(OSTYPE),linux)
	ifneq (,$(findstring 86_64, $(shell uname -m)))
		ARCH = amd64
	endif
	ifneq (,$(findstring ia64, $(shell uname -m)))
		ARCH = ia64
	endif
	ifneq (,$(findstring aarch64, $(shell uname -m)))
		ARCH = aarch64
	endif
endif
ifeq ($(OSTYPE),darwin)
	ifneq (,$(findstring x86_64, $(shell uname -m)))
		ARCH = x86_64
	endif
	ifneq (,$(findstring arm64, $(shell uname -m)))
		ARCH = arm64
	endif
endif
# For Windows, we decide whether to build in 64-bit mode based on
# whether java is 32/64-bit (since these need to match)
ifeq ($(OSTYPE),cygwin)
    JAVA_VERSION_STRING = $(shell java -version 2>&1)
    ifneq (,$(findstring 64-bit, $(JAVA_VERSION_STRING)))
        ARCH = x86_64
    endif
    ifneq (,$(findstring 64-Bit, $(JAVA_VERSION_STRING)))
        ARCH = x86_64
    endif
endif

##################
# Pre-processing #
##################

# Convert scripts needed for building, in case extracted with CRLF endings
PRE_PROC := $(shell if [ "$(OSTYPE)" = "cygwin" ]; then \
  dos2unix "$(PRISM_SRC_DIR)"/scripts/*.sh; \
  dos2unix "$(CUDD_DIR)"/setup.sh; \
fi 2> /dev/null)

########
# Java #
########

# JAVA_DIR needs to be set to the location of your Java installation.
# This makefile will try to detect this automatically based on the location of the javac command.
# If this detection does not work, or you wish to override it,
# either set the variable yourself by uncommenting and/or modifying one of the lines below
# or pass a value to make directly, e.g.: make JAVA_DIR=/usr/java

# Find javac
DETECT_JAVAC = $(shell $(PRISM_SRC_DIR)/scripts/findjavac.sh 2> /dev/null)

# Find directory containing javac
ifeq ("$(DETECT_JAVAC)","")
  JAVA_DIR =
else
  JAVA_DIR = $(shell dirname "$(DETECT_JAVAC)" | sed -e 's/\/bin//' -e 's/\/Commands//')
endif

# As a backup way of detecting JAVA_DIR, run java_home
JAVA_DIR_BACKUP = $(shell \
	if [ -f /usr/libexec/java_home ]; then /usr/libexec/java_home; \
	else echo ""; fi )

#JAVA_DIR =	/usr/java
#JAVA_DIR =	/usr/java/j2sdk1.4.2
#JAVA_DIR =	/bham/java/packages/j2sdk1.4.2
#JAVA_DIR =	/cygdrive/c/java/j2sdk1.4.2
#JAVA_DIR =	/System/Library/Frameworks/JavaVM.framework

# Now we locate the JNI header files jni.h and jni_md.h
# (in fact this is the only reason we need JAVA_DIR)
JAVA_JNI_H_DIR = $(shell \
	if [ -f "$(JAVA_DIR)"/include/jni.h ]; then echo "$(JAVA_DIR)"/include; \
	elif [ -f "$(JAVA_DIR)"/Headers/jni.h ]; then echo "$(JAVA_DIR)"/Headers; \
	elif [ -f "$(JAVA_DIR_BACKUP)"/include/jni.h ]; then echo "$(JAVA_DIR_BACKUP)"/include; \
	elif [ -f "$(JAVA_DIR_BACKUP)"/Headers/jni.h ]; then echo "$(JAVA_DIR_BACKUP)"/Headers; \
	else echo ""; fi )
JAVA_JNI_MD_H_DIR = $(shell (ls "$(JAVA_JNI_H_DIR)"/jni_md.h "$(JAVA_JNI_H_DIR)"/*/jni_md.h | head -n 1 | sed 's/\/jni_md.h//') 2>/dev/null)
JAVA_INCLUDES = -I $(JAVA_JNI_H_DIR) -I $(JAVA_JNI_MD_H_DIR)

##################
# Compilers etc. #
##################

ifeq ($(OSTYPE),cygwin)
	ifeq ($(ARCH),x86_64)
		CC = /usr/bin/x86_64-w64-mingw32-gcc
		CXX = /usr/bin/x86_64-w64-mingw32-g++
	else
		CC = /usr/bin/i686-w64-mingw32-gcc
		CXX = /usr/bin/i686-w64-mingw32-g++
	endif
	JAVACC = javacc.bat
else
	CC = gcc
	CXX = g++
	JAVACC = javacc
endif
LD = $(CXX)
JAVAC = javac
JAVA = java

export CC CXX LD JAVAC JAVACC

##############
# Flags etc. #
##############

# Tell compiler to generate debug information?
# (WARNING: must not contain a % symbol!)
DEBUG = 
#DEBUG = -g

# Compiler optimisation level:
# (WARNING: must not contain a % symbol!)
OPTIMISE = -O3
#OPTIMISE =

# Compiler warnings to enable:
# (WARNING: must not contain a % symbol!)
WARNINGS = #-Wformat

# Flags for compilation/linking
# Flags to generate shared libraries
# Executable/library naming conventions
# Suffix for binary distribution directory
# (requires GNU make for conditional evaluation)

# Linux
ifeq ($(OSTYPE),linux)
	ifeq ($(ARCH),amd64)
		# Position Independent Code required on AMD64/Itanium
		CFLAGS = -m64 -fPIC -DPIC $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS =
		BINDISTSUFFIX = linux64
		BINDISTARCH = x86
	else
	ifeq ($(ARCH),ia64)
		# Position Independent Code required on AMD64/Itanium
		# Note: We omit the -m64 flag from here since it seems to be unsupported by gcc on IA64
		CFLAGS = -fPIC -DPIC $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS =
		BINDISTSUFFIX = linux64
		BINDISTARCH = ia64
	else
	ifeq ($(ARCH),aarch64)
		# Position Independent Code required on Aarch64
		CFLAGS = -fPIC -DPIC $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS =
		BINDISTSUFFIX = linux64
		BINDISTARCH = arm
	else
		CFLAGS = -m32 $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS =
		BINDISTSUFFIX = linux32
		BINDISTARCH = x86
	endif
	endif
	endif
	BIN_TARGETS=prism.linux xprism.linux
	JFLAGS := -encoding UTF8
	SHARED = -shared
	#SHARED = -G
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .so
	LIBMATH = -lm
	CLASSPATHSEP = :
endif
# Solaris
ifeq ($(OSTYPE),solaris)
	CFLAGS = -mcpu=ultrasparc $(DEBUG) $(OPTIMISE) $(WARNINGS)
	CXXFLAGS = --std=c++11 $(CFLAGS)
	LDFLAGS =
	BINDISTSUFFIX = solaris
	BINDISTARCH = solaris
	BIN_TARGETS=prism.linux xprism.linux
	JFLAGS := -encoding UTF8
	SHARED = -shared -mimpure-text
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .so
	LIBMATH = -lm
	CLASSPATHSEP = :
endif
# Cygwin
ifeq ($(OSTYPE),cygwin)
	# -DWIN32 needed for lpsolve (included e.g. in sparse) - WIN32 is only defined by mingw/gcc without --std=c++11
	# -static-libgcc -static-libstdc+ needed to build binaries that don't rely on these libraries to be present
	# -Wl,--add-stdcall-alias needed so shared libraries can be read from JNI
	# -Wl ... -lpthread ... --no-whole-archive needed to statically include pthread in case missing
	ifeq ($(ARCH),x86_64)
		CFLAGS = $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 -DWIN32 $(CFLAGS)
		LDFLAGS = -static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive
		BINDISTSUFFIX = win64
		BINDISTARCH = x86
	else
		CFLAGS = -march=i686 $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 -DWIN32 $(CFLAGS)
		LDFLAGS = -static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive
		BINDISTSUFFIX = win32
		BINDISTARCH = x86
	endif
	BIN_TARGETS=prism.cygwin xprism.linux prism.bat.win xprism.bat.win
	JFLAGS := -encoding UTF8
	SHARED = -shared
	#SHARED = -G
	EXE = .exe
	LIBPREFIX =
	LIBSUFFIX = .dll
	LIBMATH = 
	CLASSPATHSEP = ;
endif
# Darwin
ifeq ($(OSTYPE),darwin)
	ifeq ($(ARCH),x86_64)
		CFLAGS = -arch x86_64 -fPIC -DPIC $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS = -Wl,-search_paths_first
		BINDISTSUFFIX = mac64
		BINDISTARCH = x86
		BIN_TARGETS=prism.darwin xprism.linux
	else
	ifeq ($(ARCH),arm64)
		CFLAGS = -arch arm64 -fPIC -DPIC $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS = -Wl,-search_paths_first
		BINDISTSUFFIX = mac64
		BINDISTARCH = arm
		BIN_TARGETS=prism.darwin xprism.linux
	else
		CFLAGS = -arch i386 $(DEBUG) $(OPTIMISE) $(WARNINGS)
		CXXFLAGS = --std=c++11 $(CFLAGS)
		LDFLAGS = -Wl,-search_paths_first
		BINDISTSUFFIX = mac32
		BINDISTARCH = x86
		BIN_TARGETS=prism.darwin xprism.linux
	endif
	endif
	JFLAGS := -encoding UTF8
	SHARED = -dynamiclib
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .dylib
	LIBMATH = -lm
	CLASSPATHSEP = :
endif

export CFLAGS CXXFLAGS LDFLAGS JFLAGS LIBPREFIX LIBSUFFIX

##########################################
# Main part of Makefile: Compiling PRISM #
##########################################

MAKE_DIRS = dd jdd odd dv prism mtbdd sparse hybrid parser settings userinterface pepa/compiler simulator jltl2ba jltl2dstar explicit pta param strat automata common cex

EXT_PACKAGES = lpsolve55 lp_solve_5.5_java z3 yices ppl

.PHONY: clean javadoc tests release

# inhibit building in parallel (-j option)
.NOTPARALLEL:

default: all

all: cuddpackage extpackages prism

# Build our customised version of CUDD
# (let CUDD choose its own linker - problems on Cygwin if not)
# (should be no need for linker options - we use static CUDD libraries)
# (note CFLAGS/CXXFLAGS _are_ passed on, e.g. to preserve PIC consistency)
cuddpackage: checks
	@echo Making cudd ...; cd $(CUDD_DIR) && \
	LD= LDFLAGS= ./install.sh

# Use this to force a rebuild (with javacc) of the main parser
parser:
	@echo Making parser ...; \
	cd $(PRISM_SRC_DIR)/parser && \
	$(MAKE) touch && \
	$(MAKE) \
	CLASSPATHSEP="$(CLASSPATHSEP)"

# Build various external libraries needed by PRISM
extpackages: checks
	@for ext in $(EXT_PACKAGES); do \
	  echo Making $$ext ...; \
	  (cd ext/$$ext && \
	  $(MAKE) \
	  OSTYPE="$(OSTYPE)" \
	  ARCH="$(ARCH)" \
	  SHARED="$(SHARED)" \
	  LIBMATH="$(LIBMATH)" \
	  BINDISTSUFFIX="$(BINDISTSUFFIX)" \
	  BINDISTARCH="$(BINDISTARCH)" \
	  JAVA_DIR="$(JAVA_DIR)" \
	  JAVA_JNI_H_DIR="$(JAVA_JNI_H_DIR)" \
	  JAVA_JNI_MD_H_DIR="$(JAVA_JNI_MD_H_DIR)" \
	  ) || exit 1; \
	done

# Compile main PRISM code
# (we also do some preparatory checks, and build launch scripts afterwards)
prism: checks make_dirs bin_scripts

# Compile each (top-level) source directory separately
make_dirs:
	@mkdir -p bin $(PRISM_CLASSES_DIR) $(PRISM_OBJ_DIR)/dd $(PRISM_OBJ_DIR)/jdd $(PRISM_OBJ_DIR)/odd $(PRISM_OBJ_DIR)/dv $(PRISM_OBJ_DIR)/prism $(PRISM_OBJ_DIR)/mtbdd $(PRISM_OBJ_DIR)/sparse $(PRISM_OBJ_DIR)/hybrid $(PRISM_OBJ_DIR)/simulator
	@for dir in $(MAKE_DIRS); do \
	  echo Making $(PRISM_SRC_DIR)/$$dir ...; \
	  (cd $(PRISM_SRC_DIR)/$$dir && \
	  $(MAKE) \
	  CUDD_DIR="$(CUDD_DIR)" \
	  JAVA_INCLUDES="$(JAVA_INCLUDES)" \
	  JAVA_JNI_H_DIR="$(JAVA_JNI_H_DIR)" \
	  JAVA_JNI_MD_H_DIR="$(JAVA_JNI_MD_H_DIR)" \
	  SHARED="$(SHARED)" \
	  EXE="$(EXE)" \
	  LIBMATH="$(LIBMATH)" \
	  CLASSPATHSEP="$(CLASSPATHSEP)") \
	  || exit 1; \
	done
# On Windows, convert the generated JNI headers using dos2unix
# and any scripts to be run from Cywgin, in case extracted with CRLF endings.
	@if [ "$(OSTYPE)" = "cygwin" ]; then \
	  dos2unix $(PRISM_INCLUDE_DIR)/jni/*.h; \
	  find .. -name '*.sh' -exec dos2unix {} \; ; \
	  dos2unix "$(PRISM_SRC_DIR)"/bin/prism.cygwin; \
	fi

# Compile unit tests
make_tests:
	@echo Making $(PRISM_TESTS_DIR) ...; \
	(cd $(PRISM_TESTS_DIR) && \
	$(MAKE) \
	CLASSPATHSEP="$(CLASSPATHSEP)") \
	|| exit 1;

# Copy/modify the launch scripts and put in the bin directory
bin_scripts:
	@for target in $(BIN_TARGETS); do \
	  target_trunc=`echo $$target | sed 's/\.[^.]*$$//'` && \
	  echo Copying "$(PRISM_SRC_DIR)/bin/$$target -> bin/$$target_trunc" && \
	  cp $(PRISM_SRC_DIR)/bin/$$target bin/$$target_trunc; \
	done;
	@echo Copying "$(PRISM_OBJ_DIR)/prism/ngprism$(EXE) -> bin/ngprism$(EXE)" && \
	cp $(PRISM_OBJ_DIR)/prism/ngprism$(EXE) bin/ngprism$(EXE)
	@./install.sh silent

# Some checks to make sure that the main settings are valid
checks:
	@(if [ "$(OSTYPE)" != "linux" -a "$(OSTYPE)" != "solaris" -a "$(OSTYPE)" != "cygwin" -a "$(OSTYPE)" != "darwin" ]; then \
	  echo "\033[33mTo compile PRISM, the environment variable OSTYPE"; \
	  echo "must be set to one of: linux, solaris, cygwin or darwin,"; \
	  echo "depending on which operating system you are using."; \
	  echo "This is not the case on your system. Please specify"; \
	  echo "the value of OSTYPE manually to make, e.g.:"; \
	  echo; \
	  echo "  make OSTYPE=linux"; \
	  echo; \
	  echo "Alternatively, if you wish, you can set the environment"; \
	  echo "variable yourself (using setenv or export) or you"; \
	  echo "can edit the value of OSTYPE directly in the Makefile."; \
	  echo "\033[0m"; \
	  exit 1; \
	fi; \
	if [ "$(JAVA_DIR)" = "" ]; then \
	  echo "\033[33mPRISM was unable to find the directory which contains"; \
	  echo "your Java distribution. Please specify this manually to"; \
	  echo "make, as in these examples:"; \
	  echo; \
	  echo "  make JAVA_DIR=/usr/java/j2sdk1.4.2"; \
	  echo "  make JAVA_DIR=\"/cygdrive/c/Program Files/Java/jdk1.4.2\""; \
	  echo; \
	  echo "See the PRISM manual for further information."; \
	  echo; \
	  echo "Alternatively, if you wish, you can set the environment"; \
	  echo "variable yourself (using setenv or export) or you"; \
	  echo "can edit the value of JAVA_DIR directly in the Makefile."; \
	  echo "\033[0m"; \
	  exit 1; \
	fi; \
	if [ ! -d "$(JAVA_DIR)" ]; then \
	  echo "\033[33mJava directory \"$(JAVA_DIR)\" does not exist."; \
	  echo "\033[0m"; \
	  exit 1; \
	fi; \
	if [ ! -f "$(JAVA_JNI_H_DIR)"/jni.h ]; then \
	  echo "\033[33mCould not locate JNI header jni.h within \"$(JAVA_DIR)\"."; \
	  echo "You may need to set JAVA_DIR by hand. See the PRISM manual for details."; \
	  echo "\033[0m"; \
	  exit 1; \
	fi; \
	if [ ! -f "$(JAVA_JNI_MD_H_DIR)"/jni_md.h ]; then \
	  echo "\033[33mCould not locate JNI header jni_md.h within \"$(JAVA_DIR)\"."; \
	  echo "You may need to set JAVA_DIR by hand. See the PRISM manual for details."; \
	  echo "\033[0m"; \
	  exit 1; \
	fi; \
	echo "VERSION: $(VERSION)"; \
	echo "OSTYPE/ARCH: $(OSTYPE) $(ARCH)"; \
	echo "JAVA_DIR: $(JAVA_DIR)"; \
	echo "JAVA_DIR_BACKUP: $(JAVA_DIR_BACKUP)"; \
	echo "JAVAC (which): "`which $(JAVAC)`; \
	echo "JAVAC (version): "$(shell $(JAVAC) -version) \
	)

# Misc: count the number of lines of code
count_loc:
	find $(PRISM_SRC_DIR) -name '*.java' -o -name '*.cc' | xargs wc -l

###########
# Testing #
###########

# Run all unit tests
unittests: make_tests
	# Provide Regex to match our test classes. If none is given, only certain test classes are excluded by default.
	$(JAVA) -jar lib/junit-platform-console-standalone.jar -cp classes --include-classname '^(Test.*|.+[.$$]Test.*|.+Tests?[.$$].+|.*Tests?)$$' -scan-classpath --details=summary

# Run a single test case from the test suite (useful quick check that the build was ok)
test:
	bin/prism etc/tests/dtmc_pctl.prism etc/tests/dtmc_pctl.prism.props -h -test

# Run a single test case from the test suite that relies on lpsolve being properly installed
testlpsolve:
	bin/prism etc/tests/test_lpsolve_mdpmo.prism etc/tests/test_lpsolve_mdpmo.prism.props -lp -test

# Run a single test case from the test suite that relies on PPL being properly installed
testppl:
	bin/prism etc/tests/test_ppl_smg.prism etc/tests/test_ppl_smg.prism.props -test

# Run a single test case from the test suite that relies on Z3 being properly installed
testz3:
	bin/prism etc/tests/test_z3_csg.prism etc/tests/test_z3_csg.prism.props -smtsolver z3 -test
	
# Run a single test case from the test suite that relies on Yices being properly installed
testyices:
	bin/prism etc/tests/test_yices_csg.prism etc/tests/test_yices_csg.prism.props -smtsolver yices -test

# Run all tests from the test suite (in ../prism-tests and ./tests)
# Optionally, extra arguments for prism-auto are picked up via variable TESTS_ARGS
tests: testslocal
	@if [ -d ../prism-tests ]; then \
	  etc/scripts/prism-auto -t -m ../prism-tests -p bin/prism --print-failures --nailgun --ngprism bin/ngprism $(TESTS_ARGS); \
	else \
	  echo "Skipping tests"; \
	fi

# Just display the command to run the test suite on this version of PRISM
# Optionally, extra arguments for prism-auto are picked up via variable TESTS_ARGS
testsecho:
	@echo etc/scripts/prism-auto -t -m ../prism-tests -p bin/prism --print-failures --nailgun --ngprism bin/ngprism $(TESTS_ARGS)

# Run local tests (in ./tests)
# Optionally, extra arguments for prism-auto are picked up via variable TESTS_ARGS
testslocal:
	@if [ -d tests ]; then \
	  etc/scripts/prism-auto -t -m tests -p bin/prism --print-failures --nailgun --ngprism bin/ngprism $(TESTS_ARGS); \
	else \
	  echo "Skipping local tests"; \
	fi

# Run the tests from ../prism-tests (with different engine settings, picked up from ../prism-tests/all-engines.args)
# - Export tests are disabled, as there is currently no robust test mechanism for dealing with the
#   variations in the output of the different engines.
# - We run with --test-all, as failures for some engines should not abort the tests
# - We run with a timeout of 1 minute, as some engines take a long time for some properties
testsfull:
	etc/scripts/prism-auto -t -m ../prism-tests \
	--skip-export-runs --skip-duplicate-runs --test-all -a ../prism-tests/all-engines.args --timeout 1m \
	-p bin/prism --print-failures --nailgun $(TESTS_ARGS);

##########################
# Building distributions #
##########################

# Build prism.jar
binary:
	@echo "Generating JAR file ($(PRISM_LIB_DIR)/prism.jar)..."
	@jar cmf $(PRISM_SRC_DIR)/manifest.txt $(PRISM_LIB_DIR)/prism.jar -C $(PRISM_CLASSES_DIR) . -C . $(PRISM_IMAGES_DIR) $(PRISM_DTDS_DIR)

# Build prism-sources.jar
source-jar:
	@echo "Generating sources JAR file ($(PRISM_LIB_DIR)/prism-sources.jar)..."
	@find $(PRISM_SRC_DIR) -type f -name '*.java' -o -name '*.form' -o -name '*.jj' | sed -e "s/^$(PRISM_SRC_DIR)./-C $(PRISM_SRC_DIR) /" > prism-sources.txt
	@jar cf $(PRISM_LIB_DIR)/prism-sources.jar @prism-sources.txt
	@rm -f prism-sources.txt

# Download a local html copy of the manual
#PRISM_MANUAL_WEBSITE = http://prismmodelchecker.localhost/manual/
PRISM_MANUAL_WEBSITE = http://www.prismmodelchecker.org/manual/
doc: clean_doc
	(cd .. && wget -r -np -k -E -nH --no-cookies --header "Cookie: setskin=offline" --restrict-file-names=windows --reject '*action=sourceblock*' $(PRISM_MANUAL_WEBSITE) $(PRISM_MANUAL_WEBSITE)/pub/skins/offline/images/)
clean_doc:
	rm -rf ../manual

# Set up version, in particular for building releases
# Unless VERSION has been passed in (as VERSION=xxx),
# extract version number from Java code using printversion
VERSION = $(shell $(PRISM_SRC_DIR)/scripts/printversion.sh $(PRISM_SRC_DIR) 2> /dev/null)

# Default branch for building source releases
BRANCH = master

# File/dir names for distribution
SRC_DIST_DIR = prism-games-$(VERSION)-src
SRC_DIST_FILE = $(SRC_DIST_DIR).tar.gz
BIN_DIST_DIR = prism-games-$(VERSION)-$(BINDISTSUFFIX)-$(BINDISTARCH)
BIN_DIST_FILE = $(BIN_DIST_DIR).tar.gz
BIN_DIST_INST = $(BIN_DIST_DIR)-installer.exe

# Build a (binary) distribution in the release directory

release: release_config clean_all all binary release_check_version build_release clean_binary

release_config:
	$(eval JFLAGS := --release 9 $(JFLAGS))

release_check_version:
	@if [ "$(VERSION)" = "" ]; then echo "Usage: make release VERSION=4.5"; exit 1; fi

build_release:
	@if [ "$(BINDISTSUFFIX)" = "win32" -o "$(BINDISTSUFFIX)" = "win64" ]; then \
		echo Creating Windows installer in "release/$(BIN_DIST_INST)"... && \
		makensis /NOCD /DPRISM_NAME="PRISM-games $(VERSION)" /DPRISM_BUILD="prism-games-$(VERSION)" /DPRISM_BINDISTSUFFIX="$(BINDISTSUFFIX)" /DPRISM_BINDISTARCH="$(BINDISTARCH)" /DPRISM_DIR=".." $(PRISM_SRC_DIR)/nsis_script.nsi && \
		mkdir -p release && \
		mv ../"$(BIN_DIST_INST)" release; \
	else \
		echo Creating binary distribution in "release/$(BIN_DIST_FILE)"... && \
		mkdir -p release && \
		rm -rf "release/$(BIN_DIST_DIR)" && \
		rm -f "release/$(BIN_DIST_FILE)" && \
		mkdir -p "release/$(BIN_DIST_DIR)" && \
		tar cf - bin $(PRISM_LIB_DIR) install.sh etc -C .. manual prism-examples CHANGELOG.txt COPYING.txt README.md | tar xf - -C "release/$(BIN_DIST_DIR)" && \
		(cd release && tar cfz "$(BIN_DIST_FILE)" "$(BIN_DIST_DIR)") && rm -rf "release/$(BIN_DIST_DIR)"; \
	fi

release_source: release_check_version build_release_source

build_release_source:
	@echo Creating source distribution from branch "$(BRANCH)" in "release/$(SRC_DIST_FILE)"... && \
	mkdir -p release && \
	rm -rf "release/$(SRC_DIST_DIR)" && \
	rm -f "release/$(SRC_DIST_FILE)" && \
	mkdir -p "release/$(SRC_DIST_DIR)" && \
	(cd .. && git archive $(BRANCH) --format=tar) | tar xf - -C "release/$(SRC_DIST_DIR)" && \
	(cd release && tar cfz "$(SRC_DIST_FILE)" "$(SRC_DIST_DIR)") && rm -rf "release/$(SRC_DIST_DIR)"

# Build Javadoc (and put in javadoc directory)

PRISM_CLASSPATH = "$(PRISM_CLASSES_DIR)$(CLASSPATHSEP)$(PRISM_LIB_DIR)/*"

javadoc:
	@JAVADOC_DIRS=`echo $(MAKE_DIRS) | sed 's/\//./g' | sed 's/ /:/g'` && \mkdir -p javadoc; javadoc $(JFLAGS) -d javadoc -overview $(PRISM_SRC_DIR)/overview.html -sourcepath $(PRISM_SRC_DIR) -classpath $(PRISM_SRC_DIR)$(CLASSPATHSEP)$(PRISM_CLASSPATH) -subpackages $$JAVADOC_DIRS -exclude parser

###############
# Cleaning up #
###############

# Clean main PRISM build (not CUDD or external libs)
clean: checks
	@(for dir in $(MAKE_DIRS); do \
	  echo Cleaning $(PRISM_SRC_DIR)/$$dir ...; \
	  (cd $(PRISM_SRC_DIR)/$$dir && \
	  $(MAKE) -s EXE="$(EXE)" clean) \
	  || exit 1; \
	done; \
	find $(PRISM_CLASSES_DIR) -name '*.class' -exec rm {} \; ; \
	rm -f $(PRISM_LIB_DIR)/*jnilib; \
	rm -f $(PRISM_LIB_DIR)/prism.jar; \
	rm -f $(PRISM_LIB_DIR)/prism-sources.jar; \
	rm -f $(BIN_PRISM) $(BIN_XPRISM) $(BIN_PRISM_BAT) $(BIN_XPRISM_BAT) )

celan: clean

# Clean PRISM + CUDD and external libs
clean_all: checks clean_cudd clean_ext clean clean_tests

clean_cudd:
	@(cd $(CUDD_DIR) && ./clean.sh)

clean_ext:
	@(for ext in $(EXT_PACKAGES); do \
	  echo Cleaning $$ext ...; \
	  (cd ext/$$ext && \
	  $(MAKE) -s OSTYPE="$(OSTYPE)" ARCH="$(ARCH)" clean) \
	  || exit 1; \
	done )

clean_tests:
	@(cd $(PRISM_TESTS_DIR) && $(MAKE) clean)

# Remove just the prism.jar binary
clean_binary:
	@echo "Removing JAR file ($(PRISM_LIB_DIR)/prism.jar)..."
	@rm -f $(PRISM_LIB_DIR)/prism.jar

# Clear individual directories (sometimes useful)
clean_dd: checks
	@(cd $(PRISM_SRC_DIR)/dd && $(MAKE) -s EXE="$(EXE)" clean)
clean_jdd: checks
	@(cd $(PRISM_SRC_DIR)/jdd && $(MAKE) -s EXE="$(EXE)" clean)
clean_odd: checks
	@(cd $(PRISM_SRC_DIR)/odd && $(MAKE) -s EXE="$(EXE)" clean)
clean_dv: checks
	@(cd $(PRISM_SRC_DIR)/dv && $(MAKE) -s EXE="$(EXE)" clean)
clean_prism: checks
	@(cd $(PRISM_SRC_DIR)/prism && $(MAKE) -s EXE="$(EXE)" clean)
clean_mtbdd: checks
	@(cd $(PRISM_SRC_DIR)/mtbdd && $(MAKE) -s EXE="$(EXE)" clean)
clean_sparse: checks
	@(cd $(PRISM_SRC_DIR)/sparse && $(MAKE) -s EXE="$(EXE)" clean)
clean_hybrid: checks
	@(cd $(PRISM_SRC_DIR)/hybrid && $(MAKE) -s EXE="$(EXE)" clean)
clean_parser: checks
	@(cd $(PRISM_SRC_DIR)/parser && $(MAKE) -s EXE="$(EXE)" clean)
clean_userinterface: checks
	@(cd $(PRISM_SRC_DIR)/userinterface && $(MAKE) -s EXE="$(EXE)" clean)
clean_simulator: checks
	@(cd $(PRISM_SRC_DIR)/simulator && $(MAKE) -s EXE="$(EXE)" clean)
clean_jltl2ba: checks
	@(cd $(PRISM_SRC_DIR)/jltl2ba && $(MAKE) -s EXE="$(EXE)" clean)
clean_jltl2dstar: checks
	@(cd $(PRISM_SRC_DIR)/jltl2dstar && $(MAKE) -s EXE="$(EXE)" clean)
clean_explicit: checks
	@(cd $(PRISM_SRC_DIR)/explicit && $(MAKE) -s EXE="$(EXE)" clean)
clean_pta: checks
	@(cd $(PRISM_SRC_DIR)/pta && $(MAKE) -s EXE="$(EXE)" clean)
clean_param: checks
	@(cd $(PRISM_SRC_DIR)/param && $(MAKE) -s EXE="$(EXE)" clean)
clean_strat: checks
	@(cd $(PRISM_SRC_DIR)/strat && $(MAKE) -s EXE="$(EXE)" clean)

#################################################
