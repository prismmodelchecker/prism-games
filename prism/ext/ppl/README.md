This directory contains pre-compiled versions of the Parma Polyhedra Library (PPL) library
and its Java interface. The Makefile copies the correct versions to the PRISM lib directory.

---

Instructions to build the libraries from source on different OSs are below.

We build GMP from source first. This allows us to statically link GMP
to the PPL shared libraries, minimising dependencies to bundle.
Also, we need JNI friendly libraries for Cygwin.
We use a custom version of PPL that has some build fixes and additions,
including the option to statically link GMP.

---

### Linux

(assuming set up as in prism-install-ubuntu)

```
# Pre-requisites
sudo apt -y install autoconf automake libtool

# Set-up
export BUILD_DIR="$HOME/tools"
export JAVA_HOME=/usr/lib/jvm/default-java

# Build a dynamic GMP and a static GMP (incl. C++)
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-cxx --disable-static --enable-shared --prefix=$BUILD_DIR/dynamic_gmp CC=gcc ABI=64 CFLAGS='-fPIC' CPPFLAGS=-DPIC
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-cxx --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp CC=gcc ABI=64 CFLAGS='-fPIC' CPPFLAGS=-DPIC
make && make install
# make check

# Build (custom) PPL with Java interface
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/prismmodelchecker/ppl.git
cd ppl
autoreconf -i
./configure --enable-interfaces=java --with-java="$JAVA_HOME" --disable-documentation --prefix=$BUILD_DIR --with-gmp-include="$BUILD_DIR/static_gmp/include" --with-gmp-lib-static="$BUILD_DIR/static_gmp/lib/libgmpxx.a $BUILD_DIR/static_gmp/lib/libgmp.a" CXXFLAGS="-std=c++11" JAVACFLAGS="--release 9"
make
make install
cp $BUILD_DIR/lib/libppl.so.15 $BUILD_DIR/lib/ppl/libppl_java.so $BUILD_DIR/lib/ppl/ppl_java.jar $BUILD_DIR
```

Files `libppl.so.15`, `libppl_java.so` an `ppl_java.jar` are then in `$BUILD_DIR`.

### macOS

```
# Pre-requisites
brew install autoconf automake libtool

# Set-up
export BUILD_DIR="$HOME/tools/tmp"
export JAVA_HOME=/opt/homebrew/Cellar/openjdk/20.0.2/libexec/openjdk.jdk/Contents/Home

# Build a dynamic GMP and a static GMP (incl. C++)
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-cxx --disable-static --enable-shared --prefix=$BUILD_DIR/dynamic_gmp CC=gcc ABI=64 CFLAGS='-fPIC -m64' CPPFLAGS=-DPIC
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-cxx --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp CC=gcc ABI=64 CFLAGS='-fPIC -m64' CPPFLAGS=-DPIC
make && make install
# make check

# Build (custom) PPL with Java interface
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/prismmodelchecker/ppl.git
cd ppl
autoreconf -i
./configure --enable-interfaces=java --with-java="$JAVA_HOME" --disable-documentation --prefix=$BUILD_DIR --with-gmp-include="$BUILD_DIR/static_gmp/include" --with-gmp-lib-static="$BUILD_DIR/static_gmp/lib/libgmpxx.a $BUILD_DIR/static_gmp/lib/libgmp.a" CXXFLAGS="-std=c++11" JAVACFLAGS="--release 9"
make
make install
cp $BUILD_DIR/lib/libppl.15.dylib $BUILD_DIR
cp $BUILD_DIR/lib/ppl/libppl_java.jnilib $BUILD_DIR/libppl_java.dylib
cp $BUILD_DIR/lib/ppl/ppl_java.jar $BUILD_DIR
```

Files `libppl.15.dylib`, `libppl_java.dylib` an `ppl_java.jar` are then in `$BUILD_DIR`.

---

### Cygwin

(assuming set up as in prism-install-win.bat/prism-install-cygwin)

```
# Pre-requisites
/cygdrive/c/Users/Administrator/setup-x86_64 -P autoconf -P automake -P libtool -q

# Set-up
export BUILD_DIR="/tools"
export JAVA_HOME="/cygdrive/c/Program Files/Eclipse Adoptium/jdk-11.0.21.9-hotspot"
export LIBWINPTHREAD_DIR="/cygdrive/c/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin"

# Java install without spaces (PPL configure breaks)
mkdir -p "$BUILD_DIR" && ln -s "$JAVA_HOME" "$BUILD_DIR/java"
export JAVA_HOME2="$BUILD_DIR/java"

# Build a dynamic GMP and a static GMP (incl. C++)
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --host=x86_64-w64-mingw32 --build=i686-pc-cygwin --enable-cxx LDFLAGS="-Wl,--add-stdcall-alias" --enable-shared --disable-static --prefix=$BUILD_DIR/dynamic_gmp
make && make install
# make check
mkdir -p $BUILD_DIR && cd $BUILD_DIR && mkdir -p dynamic_gmp && mkdir -p static_gmp
rm -rf gmp-6*
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --host=x86_64-w64-mingw32 --build=i686-pc-cygwin --enable-cxx LDFLAGS="-Wl,--add-stdcall-alias" --enable-static --disable-shared --prefix=$BUILD_DIR/static_gmp
make && make install
# make check

# Build (custom) PPL with Java interface
mkdir -p $BUILD_DIR && cd $BUILD_DIR
git clone https://github.com/prismmodelchecker/ppl.git
cd ppl
autoreconf -i -f
# NB: Omit --build=i686-pc-cygwin from configure otherwise jlong test fails
./configure --host=x86_64-w64-mingw32 --enable-interfaces=java --with-java="$JAVA_HOME2" --disable-documentation --prefix=$BUILD_DIR --with-gmp-include="$BUILD_DIR/static_gmp/include" --with-gmp-lib-static="$BUILD_DIR/static_gmp/lib/libgmpxx.a $BUILD_DIR/static_gmp/lib/libgmp.a" CXXFLAGS="-std=c++11" LDFLAGS="-static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive -L$LIBWINPTHREAD_DIR" JAVACFLAGS="--release 9"
make
make install
# Manually rebuild shared libs, letting g++ take care of linking and doing static gcc/c++/lpthread properly
(cd src && x86_64-w64-mingw32-g++ -shared -static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive .libs/assertions.o .libs/Box.o .libs/checked.o .libs/Checked_Number.o .libs/Float.o .libs/fpu-ia32.o .libs/BDS_Status.o .libs/Box_Status.o .libs/Og_Status.o .libs/Concrete_Expression.o .libs/Constraint.o .libs/Constraint_System.o .libs/Congruence.o .libs/Congruence_System.o .libs/Generator_System.o .libs/Grid_Generator_System.o .libs/Generator.o .libs/Grid_Generator.o .libs/Handler.o .libs/Init.o .libs/Coefficient.o .libs/Linear_Expression.o .libs/Linear_Expression_Impl.o .libs/Linear_Expression_Interface.o .libs/Linear_Form.o .libs/Scalar_Products.o .libs/MIP_Problem.o .libs/PIP_Tree.o .libs/PIP_Problem.o .libs/Poly_Con_Relation.o .libs/Poly_Gen_Relation.o .libs/BHRZ03_Certificate.o .libs/H79_Certificate.o .libs/Grid_Certificate.o .libs/Partial_Function.o .libs/Polyhedron_nonpublic.o .libs/Polyhedron_public.o .libs/Polyhedron_chdims.o .libs/Polyhedron_widenings.o .libs/C_Polyhedron.o .libs/NNC_Polyhedron.o .libs/Grid_nonpublic.o .libs/Grid_public.o .libs/Grid_chdims.o .libs/Grid_widenings.o .libs/BD_Shape.o .libs/Octagonal_Shape.o .libs/Pointset_Powerset.o .libs/CO_Tree.o .libs/Sparse_Row.o .libs/Dense_Row.o .libs/Bit_Matrix.o .libs/Bit_Row.o .libs/Ph_Status.o .libs/Grid_Status.o .libs/Variable.o .libs/Variables_Set.o .libs/Grid_conversion.o .libs/Grid_simplify.o .libs/set_GMP_memory_alloc_funcs.o .libs/stdiobuf.o .libs/c_streambuf.o .libs/globals.o .libs/mp_std_bits.o .libs/Weight_Profiler.o .libs/version.o .libs/termination.o .libs/wrap_string.o .libs/Time.o .libs/Watchdog.o .libs/Threshold_Watcher.o $BUILD_DIR/static_gmp/lib/libgmpxx.a $BUILD_DIR/static_gmp/lib/libgmp.a -g -O2 -o .libs/libppl-15.dll -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/libppl.dll.a)
(cd interfaces/Java/jni && x86_64-w64-mingw32-g++ -shared -static-libgcc -static-libstdc++ -Wl,--add-stdcall-alias -Wl,-Bstatic,--whole-archive -lpthread -Wl,-Bdynamic,--no-whole-archive .libs/ppl_java_common.o .libs/ppl_java_globals.o .libs/ppl_java_Termination.o .libs/ppl_java_Polyhedron.o .libs/ppl_java_Grid.o .libs/ppl_java_Rational_Box.o .libs/ppl_java_BD_Shape_mpz_class.o .libs/ppl_java_BD_Shape_mpq_class.o .libs/ppl_java_Octagonal_Shape_mpz_class.o .libs/ppl_java_Octagonal_Shape_mpq_class.o .libs/ppl_java_Constraints_Product_C_Polyhedron_Grid.o .libs/ppl_java_Pointset_Powerset_C_Polyhedron.o .libs/ppl_java_Pointset_Powerset_NNC_Polyhedron.o .libs/ppl_java_Double_Box.o .libs/ppl_java_BD_Shape_double.o .libs/ppl_java_Octagonal_Shape_double.o ../../../src/.libs/libppl.dll.a $BUILD_DIR/static_gmp/lib/libgmpxx.a $BUILD_DIR/static_gmp/lib/libgmp.a -g -O2 -Wl,--kill-at -o .libs/libppl_java.dll -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/libppl_java.dll.a)
make install
cp $BUILD_DIR/bin/libppl-15.dll $BUILD_DIR
cp $BUILD_DIR/lib/ppl/libppl_java.dll $BUILD_DIR/ppl_java.dll
cp $BUILD_DIR/lib/ppl/ppl_java.jar $BUILD_DIR
```

Files `libppl-15.dll`, `ppl_java.dll` an `ppl_java.jar` are then in `$BUILD_DIR`.
