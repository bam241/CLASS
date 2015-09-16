#!/bin/bash

option=$1
case "$option" in
  "-h" | "--help")

cat <<\_ACEOF
###############################################################
############## configures and compiles CLASS V4.1 #############
###############################################################

Usage: install.sh [VAR=VALUE] [OPTION]
Defaults for the options are specified in brackets.

Configuration:
  -h, --help         display this help and exit
Optional Features:
  --disable-OMP      do not compile with OpenMP support for evolution 
                     [default: enable for gcc version >= 4.1]
  --InstallLib-path=path     Location of made CLASS's libraries [default= $PWD/lib]

Some influential environment variables:
  CXX         C++ compiler command [default=g++]
  CXXFLAGS    C++ compiler flags, e.g. -I<include dir> if
              you have headers in a nonstandard directory <include dir>
  CPPFLAGS    C++ preprocessor flags, e.g. -D<special flag>


Report bugs to <leniau@subatech.in2p3.fr>.
(special thanks to PTO)
_ACEOF

exit ;;
esac

####### set default

IsGCCSupportOMP="no"
IsOMPEnable="yes"
OMPFLAGS=

ROOTCFLAGS=
ROOTGLIBS=
ROOTLIBS=

LIBDIR=${PWD}/lib

CXX=""
CXXFLAGS=""
CPPFLAGS=""

####### evaluate options
for option
do
  case $option in
   *=?*) ac_optarg=`expr "X$option" : '[^=]*=\(.*\)'` ;;
   *=)   ac_optarg= ;;
   *)    ac_optarg=yes ;;
  esac
#  echo $ac_optarg
  case $option in
  	--InstallLib-path=*)
  		LIBDIR="$ac_optarg" ;;
	--disable-OMP)
		IsOMPEnable="disable" ;;
	 *) 
	 	echo "Unrecognized option $option"
	   	exit ;;
	esac
done

####### ROOT Support
if [ -f $ROOTSYS/bin/root-config ]
then
	echo "Checking for ROOT cern lib... yes"
	ROOTCFLAGS='$(shell ${ROOTSYS}/bin/root-config --cflags)'
	ROOTGLIBS='$(shell ${ROOTSYS}/bin/root-config --glibs)'
	ROOTLIBS='$(shell ${ROOTSYS}/bin/root-config --libs)'

	if [ "$ROOTSYS/bin/root-config --features-tmva" = "no" ]
	then
		echo "TMVA is not activated : consider rebuild ROOT activating this feature"
		exit 0
	fi
	
else
		echo "Checking for ROOT cern lib... no"
		echo "********************* ERROR *********************"
		echo "** CLASS need ROOT (cern) to work"
		echo "** Either set the ROOTSYS env variable or intall"
		echo "** the ROOT library FROM SOURCES with the same C++"
		echo "** compiler you will use for comipling CLASS"
		echo "***************************************************"
		exit 0
fi

####### OMP support
### write a testconf prog
echo "#include <omp.h>" > testconf.cxx
echo "int main(){ return 0;}" >> testconf.cxx

if [ "$IsOMPEnable" = "yes" ]
then

	g++ -o testconf -fopenmp  testconf.cxx -lgomp > script.errors 2>&1
 	
	ok=`cat script.errors|wc -l`
 	if [ $ok = 0 ]
 	then
 		IsGCCSupportOMP="yes"
	else
		IsOMPEnable="no"
		IsGCCSupportOMP="no"
 	fi
	rm -f script.errors
fi

if [ "$IsOMPEnable" = "yes" ]
then
	echo "Checking for omp.h... yes"
	echo "   You can disable the use of this library with \"--disable-OMP\" option"
	OMPFLAGS="-fopenmp -DOpenMP"
	OMPLIB=-lgomp
else 
	if [ "$IsOMPEnable" = "no" ]
	then
		echo "Checking for omp.h... no"
		echo "   OpenMP support not found."
		echo "   Either gcc is to old either install libgomp.so"
	else
		if [ "$IsOMPEnable" = "disable" ]
		then
			echo "Checking for omp.h... disable"
		fi
	fi
fi

rm -f script.errors testconf testconf.cxx


#######################################################
### Preprocessor flags & Makefile variables 
#######################################################
mkdir -p config
MAKEFILE_INC="config/Makefile.config"
rm config/Makefile.config
echo "#####################">> $MAKEFILE_INC
echo "###### OPEN MP ######">> $MAKEFILE_INC
echo "#####################">> $MAKEFILE_INC
echo "OMPFLAGS=$OMPFLAGS" >> $MAKEFILE_INC
echo "OMPLIB=$OMPLIB" >> $MAKEFILE_INC
echo >> $MAKEFILE_INC
echo "#####################">> $MAKEFILE_INC
echo "####### ROOT ########">> $MAKEFILE_INC
echo "#####################">> $MAKEFILE_INC
echo "ROOTCFLAGS:=$ROOTCFLAGS" >> $MAKEFILE_INC
echo "ROOTGLIBS:=$ROOTGLIBS" >> $MAKEFILE_INC
echo "ROOTLIBS:=$ROOTLIBS" >> $MAKEFILE_INC
echo >> $MAKEFILE_INC
echo "#####################">> $MAKEFILE_INC
echo "##### COMPILER ######">> $MAKEFILE_INC
echo "#####################">> $MAKEFILE_INC

if [ "$CXX" = "" ]
then
	echo "CXX= g++" >> $MAKEFILE_INC
else
	echo "CXX=$CXX" >> $MAKEFILE_INC
fi

if [ "$CXXFLAGS" = "" ]
then
	echo "CXXFLAGS=-O2 -g -fPIC -Wall -Wno-unused">> $MAKEFILE_INC
else
	echo "CXXFLAGS=$CXXFLAGS" >> $MAKEFILE_INC
fi
if [ "$CPPFLAGS" = "" ]
then
	echo "CPPFLAGS= \$(OMPFLAGS) " >> $MAKEFILE_INC
else
	echo "CPPFLAGS=$CPPFLAGS" >> $MAKEFILE_INC
fi

echo "####Installation folder of librairies##" >>$MAKEFILE_INC
echo "LIBDIR=$LIBDIR" >> $MAKEFILE_INC
echo "Building Librairies Folder @ $LIBDIR"
mkdir -p $LIBDIR

#######################################################
### Include flags 
#######################################################
INCLUDE_INC="config/config.hxx"
if [ "$IsOMPEnable" = "yes" ]
then
	echo "#include <omp.h>" > $INCLUDE_INC
else
	echo "#define omp_get_thread_num() 0" > $INCLUDE_INC
fi

#######################################################
### compile libraries 
#######################################################
echo "####################################################"
echo "######### Compilation of CLASS libraries ###########"
echo "####################################################"
cd  source/src ; make clean ; make -j 4 ; make install ; cd ../..
echo "####################################################"
echo "########## Compilation Done  #######################"
echo "####################################################"
echo "MURE libraries installed in"
echo  "----> $LIBDIR"
echo "####################################################"
echo "######### Compilation of CLASSGUI binary ###########"
echo "####################################################"
cd gui ; make clean ; make -j 4 ; cd ../
echo "####################################################"
echo "########## Compilation Done  #######################"
echo "####################################################"
echo
#######################################################
### set the pathes of DECAY data base
#######################################################
echo
echo
echo "####################################################"
echo "########## SET DECAY DATA BASES PATHES #############"
echo "####################################################"
cd DATA_BASES/DECAY/ALL/
sed -e "s%PATHTOBASE%`pwd`%" .Decay.tmp > Decay.idx
echo "-> Done"
cd -

#######################################################
### set the environement variables
#######################################################
echo
echo
echo "####################################################"
echo "########## ENVIRONEMENT VARIABLES ##################"
echo "####################################################"

MYDefaultSHELL=$(finger ${LOGNAME} | grep Shell: | awk '{print $4}')
SHELLPreference=".$(echo "$MYDefaultSHELL" | awk -F "/bin/" '{print $2}')rc"
echo "-> Your default shell is : $MYDefaultSHELL"
echo "-> Your $SHELLPreference will be edited if CLASS_PATH CLASS_include" 
echo "   and CLASS_lib aren't defined yet "
echo
echo "CHECKING LOADED ENVIRONEMENT VARIABLES "
echo

CLASS_PATH_To_Set=""
CLASS_include_To_Set=""
CLASS_lib_To_Set=""

if [ -z "$CLASS_PATH" ]; then
	echo "Not found in your loaded $SHELLPreference."
	echo "Setting variables ..."
	echo "PRESS ENTER IF DEFAULT IS CORRECT"
	read -p "====>ENTER THE PATH TO THE CLASS root folder (defalut ${PWD}) " CLASS_PATH_To_Set
	[ -z "${CLASS_PATH_To_Set}" ] && CLASS_PATH_To_Set="${PWD}"
	echo "Path of the CLASS include folder is $CLASS_PATH_To_Set"
	echo
	read -p "====>ENTER THE PATH TO THE CLASS INCLUDE (default: $CLASS_PATH_To_Set/source/include/): " CLASS_include_To_Set
	[ -z "${CLASS_include_To_Set}" ] && CLASS_include_To_Set="${CLASS_PATH_To_Set}/source/include/"
	echo "Path of the CLASS include folder is $CLASS_include_To_Set"
	echo 
	read -p "====>ENTER THE PATH WHERE CLASS LIB ARE INSTALLED (default: $LIBDIR): " CLASS_lib_To_Set
	[ -z "${CLASS_lib_To_Set}" ] && CLASS_lib_To_Set="$LIBDIR"
	echo "Path of the CLASS lib folder is $CLASS_lib_To_Set"
	echo 

	EXPORT=
	QUOTE=
	EQUAL=
	if [ "$MYDefaultSHELL" = "/bin/tcsh" ] || [ "$MYDefaultSHELL" = "/bin/csh" ] ; then
		EXPORT="setenv "
		EQUAL=" "
	else
		EXPORT="export "
		EQUAL="="
		QUOTE="\""
	fi
		echo "" >>$HOME/$SHELLPreference
		echo "##################" >> $HOME/$SHELLPreference
		echo "####CLASSV4.1#####" >> $HOME/$SHELLPreference
		echo "##################" >> $HOME/$SHELLPreference

		echo "$EXPORT CLASS_PATH$EQUAL$QUOTE$CLASS_PATH_To_Set$QUOTE" >> $HOME/$SHELLPreference
		echo "$EXPORT CLASS_include$EQUAL$QUOTE$CLASS_include_To_Set$QUOTE" >> $HOME/$SHELLPreference
		echo "$EXPORT CLASS_lib$EQUAL$QUOTE$CLASS_lib_To_Set$QUOTE" >> $HOME/$SHELLPreference

		echo "Environnment variables added in $HOME/$SHELLPreference"


else
	echo "A CLASS_PATH is already defined in your $SHELLPreference: $CLASS_PATH"

	if [ -z "$CLASS_lib" ]; then
		echo "Path to CLASS libraries is not set "
		echo "delete the CLASS_PATH set in your bashrc or tcshrc, source bashrc(tcshrc), and rerun this script"
	else
		echo "CLASS_lib is: $CLASS_lib"

	fi

	if [ -z "$CLASS_include" ]; then
		echo "Path to CLASS includes is not set "
		echo "delete the CLASS_PATH set in your bashrc or tcshrc, source bashrc(tcshrc), and rerun this script"

	else
		echo "CLASS_include is : $CLASS_include"
	fi

fi	
echo "LOADING $HOME/$SHELLPreference ... done"
echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo " Congratulations you are now able to compile your first     "
echo "               CLASS .cxx input.                            "
echo " Please read $CLASS_PATH_To_Set/documentation/Manual/USEGUIDE.pdf"                                     
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

$MYDefaultSHELL




