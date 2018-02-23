#! /bin/bash

#_____________________________________________________________________________
# usage
# call if option -h, -help or --help, display usage and quit
function usage ()
{
    cat <<\MANUAL_EOF
###############################################################
############# configures and compiles CLASS V4.1 ##############
###############################################################

Usage: install.sh [VAR=VALUE] [OPTION]
Defaults for the options are specified in brackets.

Configuration:
  -h, --help         display this help and exit
Optional Features:
  --disable-OMP      do not compile with OpenMP support for evolution 
                     [default: enable for gcc version >= 4.1]
  --InstallLib-path=path     Install location of CLASS's libraries [default= $PWD/lib]
  --InstallGui-path=path     Install location of the GUI binary [default= $PWD/gui/bin]

Some influential environment variables:
  CXX         C++ compiler command [default=g++]
  CXXFLAGS    C++ compiler flags, e.g. -I<include dir> if
              you have headers in a nonstandard directory <include dir>
  CPPFLAGS    C++ preprocessor flags, e.g. -D<special flag>


Report bugs to <nicolas.thiolliere@subatech.in2p3.fr>.
(special thanks to PTO)
MANUAL_EOF

exit 418
}

function build ()
{
    if [ ! -d "bld" ]; then
        mkdir bld
    fi
    cd bld
    cmake ..
    make -j ${J}
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
}


##############################################################################
### calls of all functions
##############################################################################

J="1"
BUILD=false
CLEAN=false

# loop on all arguments
for arg in "$@"; do
    case $arg in
        -h|-help|--help )
            usage ;;
        --build|-build|build )
            BUILD=true ;;
        --clean|-clean|clean )
            CLEAN=true ;;
        --clean-build )
            CLEAN=true;BUILD=true ;;
        -j* )
            J="${arg//-j=/}" ;;
    esac
done


if [ "${CLEAN}" = true ]; then
    clean
fi

if [ "${BUILD}" = true ]; then
    build
fi


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

