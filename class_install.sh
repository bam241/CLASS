#! /bin/bash

### CONFIGURATION ############################################################
ROOT_v="5.34/32"

OMPenable="yes"

libDir="lib"
guiDir="gui"

# flags
ROOTCFLAGS='$(shell ${ROOTSYS}/bin/root-config --cflags)'
ROOTGLIBS='$(shell ${ROOTSYS}/bin/root-config --glibs)'
ROOTLIBS='$(shell ${ROOTSYS}/bin/root-config --libs)'

OMPFLAGS="-fopenmp -DOpenMP"
OMPLIB="-lgomp"

# compilator and options
if [[ -z "$CXX" ]]; then
	CXX="g++ -std=c++11"
fi
if [[ -z "$CXXFLAGS" ]]; then
	CXXFLAGS="-O2 -g -fPIC -finline-functions"
fi
if [[ -z "$CPPFLAGS" ]]; then
	CPPFLAGS=$OMPFLAGS
fi

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


Report bugs to <leniau@subatech.in2p3.fr>.
(special thanks to PTO)
MANUAL_EOF

exit 1
}

# colors
# test file descriptor, if it's 1 so this is standard output, else we don't need colors
if [[ -t 1 ]]; then
	# sepecial color characters
	c_default="\033[0m";
	c_black="\033[30m";
	c_red="\033[31m";
	c_green="\033[32m";
	c_yellow="\033[33m";
	c_blue="\033[34m";
	c_magenta="\033[35m";
	c_cyan="\033[36m";
	c_lgray="\033[37m";
	c_dgray="\033[90m";
	c_lred="\033[91m";
	c_lgreen="\033[92m";
	c_lyellow="\033[93m";
	c_lblue="\033[94m";
	c_lmagenta="\033[95m";
	c_lcyan="\033[96m";
	c_white="\033[97m";
	c_bold="\033[1m";
else
	c_default="";
	c_black="";
	c_red="";
	c_green="";
	c_yellow="";
	c_blue="";
	c_magenta="";
	c_cyan="";
	c_lgray="";
	c_dgray="";
	c_lred="";
	c_lgreen="";
	c_lyellow="";
	c_lblue="";
	c_lmagenta="";
	c_lcyan="";
	c_white="";
	c_bold="";
fi

### ROOT support #############################################################
function root_version () {
	$ROOTSYS/bin/root-config --version
}
function root_TMVA () {
	features=$($ROOTSYS/bin/root-config --features)
	if [[ "$features" =~ "tmva" ]]; then
		echo "ok"
	else
		echo "ko"
	fi
}

function test_root () {
	# version
	if [[ $(root_version) == "${ROOT_v}" ]]; then
		echo -e "[ROOT]  version ${ROOT_v}         [ ${c_green}ok${c_default} ]"
	else
		echo -e "[ROOT]  version ${ROOT_v}         [${c_red}fail${c_default}]"
		echo -e "Please install ROOT-${ROOT_v} :\n\thttps://root.cern.ch/content/release-53432"
		exit
	fi

	# TMVA
	if [[ $(root_TMVA) == "ok" ]]; then
		echo -e "[ROOT]  TMVA feature            [ ${c_green}ok${c_default} ]"
	else
		echo -e "[ROOT]  TMVA feature            [${c_red}fail${c_default}]"
		exit
	fi
}

### OMP ######################################################################
function write_OMP_test () {
	echo "#include <omp.h>"
	echo "int main (int,char**) { return 0; }"
}

function test_OMP () {
	write_OMP_test > test_conf_OMP.cxx
	isOMPenable=$(g++ -fopenmp  test_conf_OMP.cxx -lgomp 2>&1 | wc -l)

	rm a.out test_conf_OMP.cxx
	
	if [[ $isOMPenable == 0 ]]; then
		echo -e "[OMP]   enable                  [ ${c_green}ok${c_default} ]"
	else
		echo -e "[OMP]   enable                  [${c_red}fail${c_default}]"
		exit
	fi
}

### Makefile.config ##########################################################
function write_makefileconfig () {
	echo -e "# autogenerate by $0
#####################
###### OPEN MP ######
#####################
OMPFLAGS=$OMPFLAGS
OMPLIB=$OMPLIB

#####################
####### ROOT ########
#####################
ROOTCFLAGS:=$ROOTCFLAGS
ROOTGLIBS:=$ROOTGLIBS
ROOTLIBS:=$ROOTLIBS

#####################
##### COMPILER ######
#####################
CXX= ${CXX}
CXXFLAGS= ${CXXFLAGS}

####Installation folder of librairies##
LIBDIR=${PWD}/$LIBDIR
####Installation folder of Gui##
Gui_bin_PATH=$Gui_bin_PATH
"
}
function make_dirs () {
	mkdir -p config
	mkdir -p $LIBDIR
	mkdir -p $Gui_bin_PATH
}

function makefileconfig () {
	make_dirs

	if [[ -d config ]]; then
		echo -e "[DIR]   build config dir        [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build config dir        [${c_red}fail${c_default}]"
		exit
	fi
	if [[ -d $LIBDIR ]]; then
		echo -e "[DIR]   build lib dir           [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build lib dir           [${c_red}fail${c_default}]"
		exit
	fi
	if [[ -d $Gui_bin_PATH ]]; then
		echo -e "[DIR]   build gui dir           [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build gui dir           [${c_red}fail${c_default}]"
		exit
	fi

	write_makefileconfig > config/Makefile.config

	if [[ -f config/Makefile.config ]]; then
		echo -e "[DIR]   write Makefile.config   [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   write Makefile.config   [${c_red}fail${c_default}]"
		exit
	fi
}

### config.hxx ###############################################################
function write_confighxx () {
	omp="$1"
	if [[ $omp =~ "no" ]]; then
		echo -e "#define omp_get_thread_num() 0"
	else
		echo -e "#include <omp.h>"
	fi
}

function confighxx () {
	omp="$1"
	write_confighxx $omp > config/config.hxx

	if [[ -f config/config.hxx ]]; then
		echo -e "[FILE]  write config.hxx        [ ${c_green}ok${c_default} ]"
	else
		echo -e "[FILE]  write config.hxx        [${c_red}fail${c_default}]"
		exit
	fi
}

### compilation CLASS ########################################################

# compile_class_lib
# clean, compile and make symbolic links in CLASS
function compile_class_lib () {
	make -C source clean   && echo -e "[CLASS] clean source directory  [ ${c_green}ok${c_default} ]" || echo -e "[CLASS] clean source directory  [${c_red}fail${c_default}]";
	make -C source -j 4    && echo -e "[CLASS] compile CLASS           [ ${c_green}ok${c_default} ]" || echo -e "[CLASS] compile CLASS           [${c_red}fail${c_default}]";
	make -C source install && echo -e "[CLASS] install CLASS           [ ${c_green}ok${c_default} ]" || echo -e "[CLASS] install CLASS           [${c_red}fail${c_default}]";
}

# compile_class_gui
# clean and compile CLASSGui
function compile_class_gui () {
	make -C gui clean && echo -e "[CLASS] clean gui directory     [ ${c_green}ok${c_default} ]" || echo -e "[CLASS] clean gui directory     [${c_red}fail${c_default}]";
	make -C gui -j 4  && echo -e "[CLASS] compile GUI             [ ${c_green}ok${c_default} ]" || echo -e "[CLASS] compile GUI             [${c_red}fail${c_default}]";
}

### decay data bases #########################################################
function set_decaydata () {
	decay_dir="DATA_BASES/DECAY/ALL"
	decay_databases_path="$(pwd)/${decay_dir}"
	sed -e "s%PATHTOBASE%${decay_databases_path}%" ${decay_dir}/.Decay.tmp > ${decay_dir}/Decay.idx && echo -e "[DB]    decay databases pathes  [ ${c_green}ok${c_default} ]" || echo -e "[DB]    decay databases pathes  [${c_red}fail${c_default}]"
}

### environement variables ###################################################
function write_shellrc () {
	if [[ $SHELL =~ "csh" ]]; then
		export_="setenv"
		equal_=" "
	else
		export_="export"
		equal_="="
	fi

	echo "#! $(which $SHELL)

${export_} CLASS_PATH${equal_}${PWD}
${export_} CLASS_include${equal_}\${CLASS_PATH}/source/include
${export_} CLASS_lib${equal_}\${CLASS_PATH}/lib

${export_} PATH${equal_}\${PATH}:\${CLASS_PATH}/gui/bin
${export_} LD_LIBRARY_PATH${equal_}\${LD_LIBRARY_PATH}:\${CLASS_lib}
"
}

function test_write_shellrc () {
	write_shellrc > class_env.sh
	[[ -f class_env.sh ]] && echo -e "[ENV]   write source file       [ ${c_green}ok${c_default} ]" || echo -e "[ENV]   write source file       [${c_red}fail${c_default}]";
}

### final information ########################################################
function display_final_info () {
	echo -e "${c_bold}
Congratulations you are now able to compile your first CLASS .cxx input.
Please read ${PWD}/documentation/Manual/USEGUIDE.pdf
And source ${PWD}/class_env.sh in your .$(echo ${SHELL} | awk -F / '{ print $NF }')rc to finalize installation.${c_default}
"
}

### procedure of installation ################################################
function CLASS_install () {
	omp="$1"
	test_root

	if [[ $omp =~ "no" ]]; then
		echo -e "[OMP]   checking for omp.h      [ ${c_green}no${c_default} ]"
	else
		echo -e "[OMP]   checking for omp.h      [ ${c_green}ok${c_default} ]"
		test_OMP
	fi

	makefileconfig
	confighxx $omp

	compile_class_lib
	compile_class_gui
	
	set_decaydata
	
	test_write_shellrc
	
	display_final_info
}

# loop on all arguments
for arg in "$@"; do
	case $arg in
		-h|-help|--help )
			usage ;;
		--disable-OMP )
			OMPenable="no" ;;
		--InstallLib-path=* )
			libDir="${arg//--InstallLib-path=/}" ;;
		--InstallGui-path=* )
			guiDir="${arg//--InstallGui-path=/}" ;;
	esac
done

LIBDIR=$libDir
Gui_bin_PATH=$guiDir

CLASS_install $OMPenable
