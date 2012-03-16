# Directory containing libraries
LIBDIR          = ../../lib
# Base directory containing includes for libmctal and libvalerr
EXTERNINC	= ../external
# Directory containing includes for MURE
LOCALINC	= ../include

ROOTCFLAGS    	= `root-config --cflags`
ROOTLIBS	= `root-config --libs`
ROOTGLIBS	= `root-config --glibs`


######### nothing to change from here #########
INCLUDES	= $(LOCALINC)/*.hxx
LIBNAME 	= CLASSpkg
OBJS		= CLASS.o Reactor.o TreatmentFactory.o IsotopicVector.o ZAI.o EvolutionDataBase.o EvolutiveProduct.o 
		  		  
		  
CXX           = g++
CXXFLAGS      = -O2 -g -Wall  -fPIC -fopenmp ${MCNP_VERSION} -I$(LOCALINC) $(ROOTCFLAGS) -I$(EXTERNINC)
LD            = g++
LDFLAGS       =  -g -Wall -fPIC -L$(LIBDIR) -lvalerr -lmctal -shared 

all:		$(OBJS)
		$(LD) $(LDFLAGS) $(OBJS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"

clean:
		@rm -vf $(OBJS) *~ core

.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx
