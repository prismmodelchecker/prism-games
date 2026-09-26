#!/bin/sh
# Build the SoPlex JNI wrapper (soplexj) against a SoPlex build.
#   SOPLEX_DIR=/path/to/soplex [SOPLEX_BUILD=$SOPLEX_DIR/build] [JAVA_HOME=...] sh build.sh
# Produces: soplexj.jar and libsoplexj.so (Linux) / libsoplexj.dylib (macOS), with SoPlex linked statically.
set -e
: "${SOPLEX_DIR:?set SOPLEX_DIR to the SoPlex source directory}"
SOPLEX_BUILD=${SOPLEX_BUILD:-$SOPLEX_DIR/build}
JAVA_HOME=${JAVA_HOME:-$(dirname "$(dirname "$(readlink -f "$(command -v javac)" 2>/dev/null || command -v javac)")")}
case "$(uname)" in
  Darwin) LIB=libsoplexj.dylib; JNI_MD=darwin; SHARED="-dynamiclib -install_name @rpath/libsoplexj.dylib" ;;
  *)      LIB=libsoplexj.so;    JNI_MD=linux;  SHARED="-shared -Wl,-soname,libsoplexj.so" ;;
esac
# Libraries SoPlex was built with (GMP, MPFR, zlib) are detected from its build directory
# (soplex/config.h and CMakeCache.txt); EXTRA_LIBS / EXTRA_INC can add or override flags.
EXTRA_LIBS=${EXTRA_LIBS:-}
EXTRA_INC=${EXTRA_INC:-}
CONFIG_H="$SOPLEX_BUILD/soplex/config.h"
CACHE="$SOPLEX_BUILD/CMakeCache.txt"
cache_get() { [ -f "$CACHE" ] && sed -n "s/^$1:[A-Z]*=//p" "$CACHE" | head -1; }
with() { [ -f "$CONFIG_H" ] && grep -q "^#define SOPLEX_WITH_$1" "$CONFIG_H"; }
AUTO_LIBS=""
AUTO_INC=""
if with GMP; then
  for v in GMP_LIBRARY GMPXX_LIBRARY; do l=$(cache_get $v); [ -n "$l" ] && [ "$l" != "${l%-NOTFOUND}" ] || AUTO_LIBS="$AUTO_LIBS $l"; done
  i=$(cache_get GMP_INCLUDE_DIRS); [ -n "$i" ] && AUTO_INC="$AUTO_INC -I$i"
  [ -z "$(cache_get GMP_LIBRARY)" ] && AUTO_LIBS="$AUTO_LIBS -lgmpxx -lgmp"
fi
if with MPFR; then
  l=$(cache_get MPFR_LIBRARY); if [ -n "$l" ]; then AUTO_LIBS="$AUTO_LIBS $l"; else AUTO_LIBS="$AUTO_LIBS -lmpfr"; fi
  i=$(cache_get MPFR_INCLUDE_DIRS); [ -n "$i" ] && AUTO_INC="$AUTO_INC -I$i"
fi
if with BOOST; then
  b=$(cache_get Boost_INCLUDE_DIR); [ -n "$b" ] && AUTO_INC="$AUTO_INC -I$b"
fi
echo "Detected from SoPlex build: libs [$AUTO_LIBS ] includes [$AUTO_INC ]"
EXTRA_LIBS="$AUTO_LIBS $EXTRA_LIBS"
EXTRA_INC="$AUTO_INC $EXTRA_INC"
mkdir -p build/classes
javac --release 11 -d build/classes src/java/soplex/SoPlex.java
jar cf soplexj.jar -C build/classes .
${CXX:-g++} -O2 -fPIC -std=c++14 $EXTRA_INC -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/$JNI_MD" \
  -I"$SOPLEX_DIR/src" -I"$SOPLEX_BUILD" -c src/c/soplexj.cpp -o build/soplexj.o
${CXX:-g++} $SHARED -o $LIB build/soplexj.o "$SOPLEX_BUILD/lib/libsoplex-pic.a" $EXTRA_LIBS -lz
echo "Built soplexj.jar and $LIB"
