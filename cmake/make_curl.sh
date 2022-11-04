ROOTDIR=${PWD}
CROSS_COMPILE="x86_64-w64-mingw32"

###############
### OPENSSL ###
###############
git clone https://github.com/openssl/openssl.git
mkdir -p build/openssl
cd build/openssl
${ROOTDIR}/openssl/Configure --cross-compile-prefix=${CROSS_COMPILE}- mingw64
make -j6

# move built files to new location
mkdir -p ${ROOTDIR}/include
mkdir -p ${ROOTDIR}/lib
cp -r ${ROOTDIR}/openssl/include/* ${ROOTDIR}/include
cp -r include/* ${ROOTDIR}/include
cp -r *.a ${ROOTDIR}/lib
cp -r *.dll ${ROOTDIR}/lib

#############
### CURL ####
#############
git clone https://github.com/curl/curl.git
autoreconf -fi curl

mkdir -p build/curl
export CPPFLAGS="-I${ROOTDIR}/include"
export LDFLAGS="-L${ROOTDIR}/lib"
export AR=${CROSS_COMPILE}-ar
export AS=${CROSS_COMPILE}-as
export LD=${CROSS_COMPILE}-ld
export RANLIB=${CROSS_COMPILE}-ranlib
export CC=${CROSS_COMPILE}-gcc
export NM=${CROSS_COMPILE}-nm
export LIBS="-lssl -lcrypto"
cd build/curl
${ROOTDIR}/curl/configure --prefix=${ROOTDIR}/build --target=${CROSS_COMPILE} --host=${CROSS_COMPILE} --with-openssl
make -j6

cp -r ${ROOTDIR}/curl/include/* ${ROOTDIR}/include
cp -r lib/.libs/*.a ${ROOTDIR}/lib
cp -r lib/.libs/*.dll ${ROOTDIR}/lib
