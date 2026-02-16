################################################
#  NB: This Makefile is designed to be called  #
#      from the main PRISM Makefile. It won't  #
#      work on its own because it needs        #
#      various options to be passed in         #
################################################

default: all

all: checks libfiles rpath

# Files to copy
NATIVE_LIBS =
ifdef OSTYPE
	NATIVE_LIBS = $(wildcard $(OSTYPE)/$(ARCH)/*)
endif
NATIVE_LIBS_COPY = $(NATIVE_LIBS:$(OSTYPE)/$(ARCH)/%=../../$(PRISM_LIB_DIR)/%)
COMMON_LIBS = $(wildcard common/*)
COMMON_LIBS_COPY = $(COMMON_LIBS:common/%=../../$(PRISM_LIB_DIR)/%)

# Try and prevent accidental makes (i.e. called manually, not from top-level Makefile)
checks:
	@if [ "$(LIBSUFFIX)" = "" ]; then \
	  (echo "Error: This Makefile is designed to be called from the main PRISM Makefile"; exit 1) \
	fi; 

# Copy files/links if they exist and are newer
# If a file is in both an OS-specific directory and common, the former takes precedence

libfiles: $(NATIVE_LIBS_COPY) $(COMMON_LIBS_COPY)
	@if [ "$(NATIVE_LIBS_COPY)" = "" ]; then \
	  echo "Missing pre-built native libraries for $(OSTYPE) $(ARCH)"; \
	fi

../../$(PRISM_LIB_DIR)/%: $(OSTYPE)/$(ARCH)/%
	@if [ -f "$<" ]; then \
	  echo Copying $< $@; \
	  if [ $(OSTYPE) = "cywgin" ]; then \
	    cp $< $@; \
	  else \
	    cp -R $< $@; \
	  fi \
	fi

../../$(PRISM_LIB_DIR)/%: common/%
	@if [ -f "$<" ]; then \
	  echo Copying $< $@; \
	  if [ $(OSTYPE) = "cygwin" ]; then \
	    cp $< $@; \
	  else \
	    cp -R $< $@; \
	  fi \
	fi

rpath:
	@if [ $(OSTYPE) = "darwin" ]; then \
	  install_name_tool -id "@rpath/libyices.2.dylib" ../../$(PRISM_LIB_DIR)/libyices.2.dylib; \
	  install_name_tool -id "@rpath/libyices2java.dylib" ../../$(PRISM_LIB_DIR)/libyices2java.dylib; \
	  install_name_tool -add_rpath "@loader_path/." ../../$(PRISM_LIB_DIR)/libyices2java.dylib; \
	  install_name_tool -change "/usr/local/lib/libyices.2.dylib" "@rpath/libyices.2.dylib" ../../$(PRISM_LIB_DIR)/libyices2java.dylib; \
	fi; \
	if [ $(OSTYPE) = "linux" ]; then \
	  patchelf --set-soname "libyices.so.2.6"  ../../$(PRISM_LIB_DIR)/libyices.so.2.6; \
	  patchelf --set-soname "libyices2java.so"  ../../$(PRISM_LIB_DIR)/libyices2java.so; \
	  patchelf --set-rpath '$$ORIGIN'  ../../$(PRISM_LIB_DIR)/libyices2java.so; \
	  patchelf --replace-needed "libyices.so.2.6" "libyices.so.2.6"  ../../$(PRISM_LIB_DIR)/libyices2java.so; \
	fi
# 	@if [ $(OSTYPE) = "darwin" ]; then \
# 	  otool -L ../../$(PRISM_LIB_DIR)/libyices*dylib; \
# 	fi; \
# 	if [ $(OSTYPE) = "linux" ]; then \
# 	  readelf -d ../../$(PRISM_LIB_DIR)/libyices*so; \
# 	fi

clean: checks
	# Remove installed files if present
	@if [ "$(COMMON_LIBS_COPY)" != "" ]; then \
	  echo Removing "$(COMMON_LIBS_COPY)" && rm -f "$(COMMON_LIBS_COPY)"; \
	fi
	@if [ "$(NATIVE_LIBS_COPY)" != "" ]; then \
	  echo Removing "$(NATIVE_LIBS_COPY)" && rm -f "$(NATIVE_LIBS_COPY)"; \
	fi
	# Overapproximate and remove all native libraries files for any OS
	rm -f ../../$(PRISM_LIB_DIR)/*lib*yices*
	rm -f ../../$(PRISM_LIB_DIR)/*yice*dll*

celan:	clean


#################################################
