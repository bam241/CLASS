#! /bin/bash

#_____________________________________________________________________________
# usage
# call if option -h, -help or --help, display usage and quit
function usage ()
{
echo "--------------------------------------------"
echo "--- CLASS Installation options "
echo "--------------------------------------------"
echo ""
echo "- ./install.sh [options]"
echo ""
echo "--- options = -h|-help|--help"
echo "-> Print help! "
echo ""
echo "--- options = --build|-build|build"
echo "-> Build following packages : CLASSpkg_root CLASSpkg CLASSGui "
echo ""
echo "--- options = --gtest|-gtest|gtest"
echo "-> Run google tests"
echo ""
echo "--- options = --clean|-clean|clean )"
echo "-> Clean the repo"
echo ""
echo "--- options = --clean-build"
echo "-> Clean first and then build"
echo ""
echo "--- options = --build-gtest"
echo "-> Build and run Google test"
echo ""
echo "--------------------------------------------"
exit 418
}

J="1"
BUILD=false
CLEAN=false
GTEST=false

function build ()
{
    if [ ! -d "bld" ]; then
        mkdir bld
    fi
    cd bld

        cmake ..
    if [ "$GTEST" == "false" ]; then
        make -j ${J} CLASSpkg_root CLASSpkg CLASSGui
        exit 0
    else
        make -j ${J}
        ../bin/RunTest
        exit 0
    fi
}

function clean ()
{
    if [ -d "bld" ]; then
        cd bld
        make clean
        cd ..
        rm -rf bld
    fi
    rm -rf source/src/*Dict.cxx
    rm -rf lib
    rm -rf bin
    
    if [ "$BUILD" = "false" ]; then
        exit 0
    fi

}

##############################################################################
### calls of all functions
##############################################################################


# Help if no argument
if [ -z "$*" ]; then usage; fi

# loop on all arguments
for arg in "$@"; do
    case $arg in
        -h|-help|--help )
            usage ;;
        --build|-build|build )
            BUILD=true ;;
        --gtest|-gtest|gtest )
            BUILD=true;GTEST=true ;;
        --clean|-clean|clean )
            CLEAN=true ;;
        --clean-build )
            CLEAN=true;BUILD=true ;;
        --build-gtest )
            BUILD=true;GTEST=true ;;
        -j* )
            J="${arg//-j=/}" ;;
        * )
            usage ;;
    esac
done


# Test is gtest is already there or if we have an internet connection
if [ "${GTEST}" = true ]; then
    echo -e "GET http://google.com HTTP/1.0\n\n" | nc google.com 80 > /dev/null 2>&1
    if [ ! $? -eq 0 ] && [ ! -d "bld/GTest/gtest/src/gtest" ]; then
        echo "An internet connection is required to compile the test"
        exit 1
    fi
fi


if [ "${CLEAN}" = true ]; then
    clean
fi

if [ "${BUILD}" = true ]; then
    build
fi


usage
exit 0

EOF
##############################################################################

%DOC

# Exit status
In `bash`, a correct exectution of a script return an exit status equals to `0` (that why a *C*/*C++* program end by `return 0;`). There is no standard for other exit status, so I (Josselin Massot) use the list of status code in HTTP for exit status.

* `404` : *Not Found*, a test find a file which doesn't exist.
* `418` : *I'm a teapot*, if only return help (I don't find a correct exit status for this).
* `501` : *Not Implemented*, test a feature which doesn't implemented on OS.
* `505` : *HTTP Version not supported*, the version of a library is not the one is expected.
* `507` : *Insufficient storage*, an error on making a directory, maybe because of the insufficient storage.

