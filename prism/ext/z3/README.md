This directory contains pre-compiled versions of Z3 and its Java interface.
The Makefile copies the correct versions to the PRISM lib directory.

---

For macOS/Windows, we take files directly from the Z3 binary releases at

https://github.com/Z3Prover/z3/releases

For Linux, we build from source (instructions below)

* libz3.{so,dylib,dll}
* libz3java.{so,dylib,dll}
* com.microsoft.z3.jar

---

### Linux

(assuming set up as in prism-install-ubuntu)

```
# Set-up
export BUILD_DIR="$HOME/tools"

# Build Z3 (with Java wrapper)
mkdir -p $BUILD_DIR && cd $BUILD_DIR
wget https://github.com/Z3Prover/z3/archive/z3-4.12.4.tar.gz
tar xvfz z3-4.12.4.tar.gz
cd z3-z3-4.12.4
python3 scripts/mk_make.py --java
cd build
make
#sudo make install
cp libz3.so libz3java.so com.microsoft.z3.jar $BUILD_DIR
```

Files `libz3.so`, `libz3java.so` and `com.microsoft.z3.jar` are then in `$BUILD_DIR`.
