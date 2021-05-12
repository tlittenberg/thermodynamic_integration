

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./install.sh /path/to/install/destination"
    exit 1
fi

INSTALL_PREFIX=$1

rm -rf build
mkdir -p build
pushd build/
cmake .. \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=true

cmake --build . -- VERBOSE=1
cmake --build . --target install
popd
