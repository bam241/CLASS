#! /bin/bash

### CONFIGURATION ############################################################
# ROOT_v : require version of ROOT
ROOT_v="5.34/32" # useless variable...
ROOT_v_low="5.34"
ROOT_v_up="5.99"

# flags
ROOTCFLAGS='$(shell ${ROOTSYS}/bin/root-config --cflags)'
ROOTGLIBS='$(shell ${ROOTSYS}/bin/root-config --glibs)'
ROOTLIBS='$(shell ${ROOTSYS}/bin/root-config --libs)'

OMPFLAGS="-fopenmp -DOpenMP"
OMPLIB="-lgomp"

# default value, could be change by argument (see usage for more information)
OMPenable="yes"
libDir="lib"
guiDir="bin"


#_____________________________________________________________________________
# compilator and options
if [[ -z "$CXX" ]]; then
	CXX="g++ -std=c++11" # default compiler is g++ with C++11
fi
if [[ -z "$CXXFLAGS" ]]; then
	CXXFLAGS="-O2 -g -fPIC -finline-functions" # some optimization flags for GCC, -fPIC and -finline-functions are not avalable for clang
fi
if [[ -z "$CPPFLAGS" ]]; then
	CPPFLAGS=$OMPFLAGS
fi

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


Report bugs to <leniau@subatech.in2p3.fr>.
(special thanks to PTO)
MANUAL_EOF

exit 418
}

# colors
# test file descriptor, if it's 1 so this is standard output, else we don't need colors
	# explication : si le script est exéctué normalement, le `file descriptor`
	# de sortie est la sortie standard, donc le `file descriptor` 1, sinon il
	# y a redirection du flux donc il ne faut pas afficher les caractères de
	# changement de couleur. (NB: si la redirection se fait dans le flux
	# d'erreur, il n'y a pas de couleur non plus).
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
	c_default=""; c_black=""; c_red=""; c_green=""; c_yellow=""; c_blue=""; c_magenta=""; c_cyan=""; c_lgray=""; c_dgray=""; c_lred=""; c_lgreen=""; c_lyellow=""; c_lblue=""; c_lmagenta=""; c_lcyan=""; c_white=""; c_bold="";
fi

### ROOT support #############################################################
# root_version
# return the version of ROOT without the patch release (only major.minor)
function root_version () {
	$ROOTSYS/bin/root-config --version | cut -d '/' -f 1
}
#_____________________________________________________________________________
# root_TMVA
# return "ok" if the version of ROOT implement the TMVA feature, "ko" else
function root_TMVA () {
	local features=$($ROOTSYS/bin/root-config --features)
	if [[ "$features" =~ "tmva" ]]; then
		echo "ok"
	else
		echo "ko"
	fi
}

#_____________________________________________________________________________
# test_root
# test version and TMVA feature of ROOT
function test_root () {
	# test version between ROOT_v_low and ROOT_v_up (remove . and patch release)
	if [[ "$(echo ${ROOT_v_low} | tr -d '.')" -le "$(echo $(root_version) | tr -d '.')" ]] && [[ "$(echo $(root_version) | tr -d '.')" -le "$(echo ${ROOT_v_up} | tr -d '.')" ]]; then
		echo -e "[ROOT]  version between ${ROOT_v_low} and ${ROOT_v_up}      [ ${c_green}ok${c_default} ]"
	else
		echo -e "[ROOT]  version between ${ROOT_v_low} and ${ROOT_v_up}      [${c_red}fail${c_default}]"
		echo -e "Please install ROOT-${ROOT_v} :\n\thttps://root.cern.ch/content/release-53432"
		exit 505
	fi

	# TMVA
	if [[ $(root_TMVA) == "ok" ]]; then
		echo -e "[ROOT]  TMVA feature                       [ ${c_green}ok${c_default} ]"
	else
		echo -e "[ROOT]  TMVA feature                       [${c_red}fail${c_default}]"
		exit 501
	fi
}

### OMP ######################################################################
# write_OMP_test
# write a small file to test OpenMP feature
function write_OMP_test () {
	echo "#include <omp.h>"
	echo "int main (int,char**) { return 0; }"
}

#_____________________________________________________________________________
# test_OMP
# if OpenMP is enable, test the OMP feature by compiling the small test
function test_OMP () {
	write_OMP_test > test_conf_OMP.cxx
	local isOMPenable=$(g++ -fopenmp  test_conf_OMP.cxx -lgomp  2>&1 | wc -l | awk '{ print $1 }')

	rm a.out test_conf_OMP.cxx 
	
	if [[ $isOMPenable == 0 ]]; then
		echo -e "[OMP]   enable                             [ ${c_green}ok${c_default} ]"
	else
		echo -e "[OMP]   enable                             [${c_red}fail${c_default}]"
		exit 501
	fi
}

### Makefile.config ##########################################################
# write_makefileconfig
# write Makefile.config file with correct value
function write_makefileconfig () {
	now=`date +%F\ %T`
	echo -e "# autogenerate by $0 ($now)"
	echo -e "
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
#_____________________________________________________________________________
# make-dirs
# make config, lib, gui/bin directories
function make_dirs () {
	mkdir -p config
	mkdir -p $LIBDIR
	mkdir -p $Gui_bin_PATH
}

#_____________________________________________________________________________
# makefileconfig
# test directories, write Makefile.config, an test it
function makefileconfig () {
	make_dirs

	if [[ -d config ]]; then
		echo -e "[DIR]   build config dir                   [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build config dir                   [${c_red}fail${c_default}]"
		exit 507
	fi
	if [[ -d $LIBDIR ]]; then
		echo -e "[DIR]   build lib dir                      [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build lib dir                      [${c_red}fail${c_default}]"
		exit 507
	fi
	if [[ -d $Gui_bin_PATH ]]; then
		echo -e "[DIR]   build gui dir                      [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build gui dir                      [${c_red}fail${c_default}]"
		exit 507
	fi

	write_makefileconfig > config/Makefile.config

	if [[ -f config/Makefile.config ]]; then
		echo -e "[DIR]   write Makefile.config              [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   write Makefile.config              [${c_red}fail${c_default}]"
		exit 404
	fi
}

### config.hxx ###############################################################
# write_confighxx
# wrtie config.hxx file
function write_confighxx () {
	local omp="$1"
	if [[ $omp =~ "no" ]]; then
		echo -e "#define omp_get_thread_num() 0"
	else
		echo -e "#include <omp.h>"
	fi
}

#_____________________________________________________________________________
# confighxx
# write and test config/config.hxx file
function confighxx () {
	local omp="$1"
	write_confighxx $omp > config/config.hxx

	if [[ -f config/config.hxx ]]; then
		echo -e "[FILE]  write config.hxx                   [ ${c_green}ok${c_default} ]"
	else
		echo -e "[FILE]  write config.hxx                   [${c_red}fail${c_default}]"
		exit 404
	fi
}

### compilation CLASS ########################################################

# compile_class_lib
# clean, compile and make symbolic links in CLASS
function compile_class_lib () {
	# clean
	make -C source clean   && echo -e "[CLASS] clean source directory             [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] clean source directory             [${c_red}fail${c_default}]"; exit 418)

	# links
	make -C source external_link && echo -e "[CLASS] create external link               [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] create external link               [${c_red}fail${c_default}]"; exit 418)
	
	# make dir for *.o files
	mkdir -p source/obj
	if [[ -d source/obj ]]; then
		echo -e "[DIR]   build obj dir                      [ ${c_green}ok${c_default} ]"
	else
		echo -e "[DIR]   build obj dir                      [${c_red}fail${c_default}]"
		exit 507
	fi
	
	# compile
	make -C source -j 4    && echo -e "[CLASS] compile CLASS                      [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] compile CLASS                                             [${c_red}fail${c_default}]"; exit 501)

	# make links (install)
	make -C source install && echo -e "[CLASS] install CLASS                      [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] install CLASS           [${c_red}fail${c_default}]"; exit 501)
}

#_____________________________________________________________________________
# compile_class_gui
# clean and compile CLASSGui
function compile_class_gui () {
	make -C gui clean && echo -e "[CLASS] clean gui directory                [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] clean gui directory     [${c_red}fail${c_default}]"; exit 501)
	make -C gui -j 4  && echo -e "[CLASS] compile GUI                        [ ${c_green}ok${c_default} ]" || (echo -e "[CLASS] compile GUI             [${c_red}fail${c_default}]"; exit 501)
}

### decay data bases #########################################################
# set_decaydata
# write DATA_BASE/DECAY/ALL/Decay.idx file with the correct path
function set_decaydata () {
	local decay_dir="DATA_BASES/DECAY/ALL"
	local decay_databases_path="$(pwd)/${decay_dir}"
	sed -e "s%PATHTOBASE%${decay_databases_path}%" ${decay_dir}/.Decay.tmp > ${decay_dir}/Decay.idx && echo -e "[DB]    decay databases pathes             [ ${c_green}ok${c_default} ]" || (echo -e "[DB]    decay databases pathes             [${c_red}fail${c_default}]"; exit 505)
}

### environement variables ###################################################
# write_shellrc
# write a small file with environment variables, the format depends on the default shell
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

#_____________________________________________________________________________
# test_write_shellrc
# write and test the class_env.sh file
function test_write_shellrc () {
	write_shellrc > class_env.sh
	if [[ -f class_env.sh ]]; then
		echo -e "[ENV]   write source file                  [ ${c_green}ok${c_default} ]";
	else
		echo -e "[ENV]   write source file                  [${c_red}fail${c_default}]";
		exit 404
	fi
}

### final information ########################################################
# display_final_info
# display final information to congragulate the user
function display_final_info () {
	echo -e "${c_bold}
Congratulations you are now able to compile your first CLASS .cxx input.
Please read ${PWD}/documentation/Manual/USEGUIDE.pdf
And source ${PWD}/class_env.sh in your \$HOME/.$(echo ${SHELL} | awk -F / '{ print $NF }')rc to finalize installation.${c_default}

echo \"source ${PWD}/class_env.sh\" >> \$HOME/.$(echo ${SHELL} | awk -F / '{ print $NF }')rc
"
}

### procedure of installation ################################################
# CLASS_install
# complet procedure of installation
function CLASS_install () {
	local omp="$1"
	test_root

	if [[ $omp =~ "no" ]]; then
		echo -e "[OMP]   checking for omp.h                 [ ${c_green}no${c_default} ]"
	else
		echo -e "[OMP]   checking for omp.h                 [ ${c_green}ok${c_default} ]"
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

##############################################################################
### calls of all functions
##############################################################################

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

# change default value with argument
LIBDIR=$libDir
Gui_bin_PATH=$guiDir

# call the complet install procedure
CLASS_install $OMPenable

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

