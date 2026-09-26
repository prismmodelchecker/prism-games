# SoPlex JNI wrapper (soplexj)

A minimal Java interface to the SoPlex LP solver (floating point), for solving the
many small LPs arising in CSG model checking (matrix games, correlated equilibria).
SoPlex has no official Java bindings; this wraps its C++ API directly.

* `src/java/soplex/SoPlex.java` – Java class (`soplex.SoPlex`): reusable LP object
  (`clear`, `addCol`, `addRow`, `optimize`, `getObjValue`, `getPrimal`, `getDual`, ...)
  plus `matrixGame(...)`, which builds and solves a matrix-game LP in one native call.
* `src/c/soplexj.cpp` – JNI implementation.
* `build.sh` – builds `soplexj.jar` and `libsoplexj.{so,dylib}` against a SoPlex build
  (SoPlex linked statically from `libsoplex-pic.a`).

Defaults are tuned for tiny LPs: SoPlex's presolver and timer are off (a time limit re-enables the timer).

## Building (macOS example, SoPlex built with GMP/MPFR/Boost from MacPorts)

```
JAVA_HOME=$(/usr/libexec/java_home) \
SOPLEX_DIR=$HOME/Desktop/Tools/soplex SOPLEX_BUILD=$HOME/Desktop/Tools/soplex/build \
CXX=clang++ sh build.sh
cp soplexj.jar libsoplexj.dylib ../../lib/
```

GMP, MPFR and Boost are picked up automatically from SoPlex's build directory (`soplex/config.h`,
`CMakeCache.txt`). Extra flags can be given with `EXTRA_INC` / `EXTRA_LIBS`.

## Licence

SoPlex is Apache 2.0; this wrapper is part of PRISM (GPL v2 or later).
