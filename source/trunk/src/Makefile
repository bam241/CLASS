# Directory containing libraries
LIBDIR          = ../../lib
# Base directory containing includes for libmctal and libvalerr
EXTERNINC	= ../external
# Directory containing includes for CLASS
LOCALINC	= ../include


ROOTCFLAGS    	= `root-config --cflags`
ROOTLIBS	= `root-config --libs`
ROOTGLIBS	= `root-config --glibs`


######### nothing to change from here #########
INCLUDES	= $(LOCALINC)/*.hxx
LIBNAME 	= CLASSpkg
OBJS		= 	CLASS.o \
			TreatmentFactory.o TreatmentFactoryDict.o\
			Reactor.o ReactorDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAI.o ZAIDict.o \
			EvolutionDataBase.o \
			EvolutiveProduct.o  
		  		  



CXX           = g++
CXXFLAGS      = -O2 -g -Wall -fopenmp -fPIC -I$(LOCALINC) $(ROOTCFLAGS) -I$(EXTERNINC)
LD            = g++
LDFLAGS       =  -g -Wall -fPIC $(ROOTLIBS) -shared 

all:		$(OBJS)
		$(LD) $(LDFLAGS) $(OBJS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"



ReactorDict.cxx:		$(LOCALINC)/Reactor.hxx
				rootcint -f $@ -c $^
 
TreatmentFactoryDict.cxx:		$(LOCALINC)/TreatmentFactory.hxx
				rootcint -f $@ -c $^

IsotopicVectorDict.cxx:		$(LOCALINC)/IsotopicVector.hxx
				rootcint -f $@ -c $^

ZAIDict.cxx:			$(LOCALINC)/ZAI.hxx
				rootcint -f $@ -c $^




clean:
		@rm -vf $(OBJS) *~ core *Dict.cxx

.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx
