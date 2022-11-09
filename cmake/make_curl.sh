ROOTDIR=${PWD}
CROSS_COMPILE="x86_64-w64-mingw32"

git clone --branch openssl-3.0.7 https://github.com/openssl/openssl.git src/openssl
git clone --branch curl-7_86_0 https://github.com/curl/curl.git src/curl
###############
### OPENSSL ###
###############
mkdir -p ${ROOTDIR}/build/openssl
cd build/openssl
${ROOTDIR}/src/openssl/Configure --prefix=/usr/${CROSS_COMPILE}/usr/local/ --cross-compile-prefix=${CROSS_COMPILE}- mingw64 #--prefix=${ROOTDIR}/openssl 
make -j6
sudo make install_sw

############
### CURL ###
############
autoreconf -fi ${ROOTDIR}/src/curl
mkdir -p ${ROOTDIR}/build/curl
cd ${ROOTDIR}/build/curl
${ROOTDIR}/src/curl/configure --prefix=/usr/${CROSS_COMPILE}/usr/local/ --target=${CROSS_COMPILE} --host=${CROSS_COMPILE} --with-openssl=/usr/${CROSS_COMPILE}/usr/local/openssl #--prefix=${ROOTDIR}/curl --exec-prefix ${ROOTDIR}/curl 
make -j6
sudo make install
