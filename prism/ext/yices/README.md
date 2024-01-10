This directory contains pre-compiled versions of Yices 2 and its Java interface.
The Makefile copies the correct versions to the PRISM lib directory.

---

Instructions to build the libraries from source on different OSs are below.

We build GMP and Yices 2 from source. This allows us to statically link GMP
to the Yices 2 Java wrapper (and to Yices), minimising dependencies to bundle.
Also, we need PIC libraries for macOS, and JNI friendly libraries for Cygwin.

For GMP and Yices 2, see also:

* [Yices 2 compilation instructions](https://github.com/SRI-CSL/yices2/blob/master/doc/COMPILING)
* [Yices 2 GMP compilation instructions](https://github.com/SRI-CSL/yices2/blob/master/doc/GMP)
* [Yices 2 Windows build code](https://github.com/SRI-CSL/yices2/blob/master/.github/workflows/windows_ci.yml)

For the Java wrapper, we don't use the provided scripts since they don't work
for Cygwin and don't reliably produce jars (and we need extra flags).

---

### Linux

(assuming set up as in prism-install-ubuntu)

```
# Pre-requisites
sudo apt -y install autoconf automake libtool gperf

# Set-up
export BUILD_DIR="$HOME/tools"
export JAVA_HOME=/usr/lib/jvm/default-java

# Build a dynamic GMP and a static GMP
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --disable-static --enable-shared --prefix=$BUILD_DIR/dynamic_gmp CC=gcc ABI=64 CFLAGS='-fPIC' CPPFLAGS=-DPIC
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp CC=gcc ABI=64 CFLAGS='-fPIC' CPPFLAGS=-DPIC
make && make install
# make check

# Build Yices 2
mkdir -p $BUILD_DIR && cd $BUILD_DIR
wget https://yices.csl.sri.com/releases/2.6.4/yices-2.6.4-src.tar.gz
tar xfz yices-2.6.4-src.tar.gz
cd yices2-Yices-2.6.4
autoconf
./configure CPPFLAGS=-I$BUILD_DIR/static_gmp/include LDFLAGS=-L$BUILD_DIR/static_gmp/lib --with-pic-gmp=$BUILD_DIR/static_gmp/lib/libgmp.a --with-pic-gmp-include-dir=$BUILD_DIR/static_gmp/include
make MODE=release static-dist
cp build/*/static_dist/lib/libyices.so.2.6.4 $BUILD_DIR
ln -s $BUILD_DIR/libyices.so.2.6.4 $BUILD_DIR/libyices.so.2.6
ln -s $BUILD_DIR/libyices.so.2.6 $BUILD_DIR/libyices.so.2
ln -s $BUILD_DIR/libyices.so.2 $BUILD_DIR/libyices.so
cp build/*/static_dist/include/* $BUILD_DIR

# Build Yices Java bindings
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/SRI-CSL/yices2_java_bindings
cd yices2_java_bindings
cd src/main/java/com/sri/yices
javac --release 9 -d ../../../../../../dist/lib -h . *.java
g++ -fPIC -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/linux" -I$BUILD_DIR -I$BUILD_DIR/static_gmp/include -c yicesJNI.cpp
g++ -shared -o libyices2java.so yicesJNI.o $BUILD_DIR/static_gmp/lib/libgmp.a -L$BUILD_DIR -lyices
cd ../../../../../..
jar cvfm yices.jar MANIFEST.txt -C dist/lib .
cp src/main/java/com/sri/yices/libyices2java.so yices.jar $BUILD_DIR
```

Files `libyices.so.2.6`, `libyices2java.so` and `yices.jar` are then in `$BUILD_DIR`.

### macOS

```
# Pre-requisites
brew install autoconf automake libtool gperf

# Set-up
export BUILD_DIR="$HOME/tools/tmp"
export JAVA_HOME=/opt/homebrew/Cellar/openjdk/20.0.2/libexec/openjdk.jdk/Contents/Home

# Build a dynamic GMP and a static GMP
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --disable-static --enable-shared --prefix=$BUILD_DIR/dynamic_gmp CC=gcc ABI=64 CFLAGS='-fPIC -m64' CPPFLAGS=-DPIC
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp CC=gcc ABI=64 CFLAGS='-fPIC -m64' CPPFLAGS=-DPIC
make && make install
# make check

# Build Yices 2
mkdir -p $BUILD_DIR && cd $BUILD_DIR
wget https://yices.csl.sri.com/releases/2.6.4/yices-2.6.4-src.tar.gz
tar xfz yices-2.6.4-src.tar.gz
cd yices2-Yices-2.6.4
autoconf
./configure CPPFLAGS=-I$BUILD_DIR/static_gmp/include LDFLAGS=-L$BUILD_DIR/static_gmp/lib --with-pic-gmp=$BUILD_DIR/static_gmp/lib/libgmp.a --with-pic-gmp-include-dir=$BUILD_DIR/static_gmp/include
make OPTION=64bits MODE=release static-dist
cp build/*/static_dist/lib/libyices.2.dylib $BUILD_DIR
ln -s $BUILD_DIR/libyices.2.dylib $BUILD_DIR/libyices.dylib
cp build/*/static_dist/include/* $BUILD_DIR

# Build Yices Java bindings
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/SRI-CSL/yices2_java_bindings
cd yices2_java_bindings
cd src/main/java/com/sri/yices
javac --release 9 -d ../../../../../../dist/lib -h . *.java
g++ -fPIC -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/darwin" -I$BUILD_DIR -I$BUILD_DIR/static_gmp/include -c yicesJNI.cpp
g++ -dynamiclib -o libyices2java.dylib yicesJNI.o $BUILD_DIR/static_gmp/lib/libgmp.a -L$BUILD_DIR -lyices
cd ../../../../../..
jar cvfm yices.jar MANIFEST.txt -C dist/lib .
cp src/main/java/com/sri/yices/libyices2java.dylib yices.jar $BUILD_DIR
```

Files `libyices.2.dylib`, `libyices2java.dylib` and `yices.jar` are then in `$BUILD_DIR`.

---

### Cygwin

(assuming set up as in prism-install-win.bat/prism-install-cygwin)

```
# Pre-requisites
/cygdrive/c/Users/Administrator/setup-x86_64 -P autoconf -P automake -P libtool -P gperf -q

# Set-up
export BUILD_DIR="/tools"
export JAVA_HOME="/cygdrive/c/Program Files/Eclipse Adoptium/jdk-11.0.21.9-hotspot"

# Build a dynamic GMP and a static GMP
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --host=x86_64-w64-mingw32 --build=i686-pc-cygwin LDFLAGS="-Wl,--add-stdcall-alias" --enable-shared --disable-static --prefix=$BUILD_DIR/dynamic_gmp
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --host=x86_64-w64-mingw32 --build=i686-pc-cygwin LDFLAGS="-Wl,--add-stdcall-alias" --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp
make && make install
# make check

# Build Yices 2
mkdir -p $BUILD_DIR && cd $BUILD_DIR
wget https://yices.csl.sri.com/releases/2.6.4/yices-2.6.4-src.tar.gz
tar xfz yices-2.6.4-src.tar.gz
cd yices2-Yices-2.6.4
autoconf
./configure --host=x86_64-w64-mingw32 CPPFLAGS=-I$BUILD_DIR/dynamic_gmp/include LDFLAGS="-L$BUILD_DIR/dynamic_gmp/lib -Wl,--add-stdcall-alias" --with-static-gmp=$BUILD_DIR/static_gmp/lib/libgmp.a --with-static-gmp-include-dir=$BUILD_DIR/static_gmp/include
export LD_LIBRARY_PATH=/usr/local/lib/:${LD_LIBRARY_PATH}
make OPTION=mingw64 MODE=release static-dist
cp build/*/static_dist/bin/libyices.dll $BUILD_DIR
cp build/*/static_dist/include/* $BUILD_DIR

# Build Yices Java bindings
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/SRI-CSL/yices2_java_bindings
cd yices2_java_bindings
cd src/main/java/com/sri/yices
javac --release 9 -d ../../../../../../dist/lib -h . *.java
x86_64-w64-mingw32-g++ -I$JAVA_HOME/include -I$JAVA_HOME/include/win32 -I$BUILD_DIR -I$BUILD_DIR/static_gmp/include -fpermissive -c yicesJNIforWindows.cpp
x86_64-w64-mingw32-g++ -shared -static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive -o yices2java.dll yicesJNIforWindows.o $BUILD_DIR/static_gmp/lib/libgmp.a -L$BUILD_DIR -lyices
cd ../../../../../..
jar cvfm yices.jar MANIFEST.txt -C dist/lib .
cp src/main/java/com/sri/yices/yices2java.dll yices.jar $BUILD_DIR
```

Files `libyices.dll`, `yices2java.dll` and `yices.jar` are then in `$BUILD_DIR`.
