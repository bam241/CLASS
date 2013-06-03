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
			Storage.o StorageDict.o\
			FabricationPlant.o FabricationPlantDict.o \
			Pool.o PoolDict.o\
			Reactor.o ReactorDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAI.o ZAIDict.o \
			DataBank.o \
			EvolutionData.o EvolutionDataDict.o \
			LogFile.o

ROOTOBJS	= 	Storage.o StorageDict.o\
			Pool.o PoolDict.o\
			Reactor.o ReactorDict.o \
			FabricationPlant.o FabricationPlantDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAI.o ZAIDict.o \
			DataBank.o \
			EvolutionData.o EvolutionDataDict.o \
			LogFile.o


CXX           = g++
CXXFLAGS      = -O2 -g -Wall -fopenmp -fPIC -I$(LOCALINC) $(ROOTCFLAGS)
LD            = g++
LDFLAGS       =  -g -Wall -fPIC $(ROOTLIBS) -shared -lgomp

all:		$(OBJS)
		$(LD) $(LDFLAGS) $(OBJS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"
		$(LD) $(LDFLAGS) $(ROOTOBJS) -o $(LIBDIR)/lib$(LIBNAME)_root.so
		@echo "lib$(LIBNAME)_root.so done"




StorageDict.cxx:		$(LOCALINC)/Storage.hxx
				rootcint -f $@ -c $^

ReactorDict.cxx:		$(LOCALINC)/Reactor.hxx
				rootcint -f $@ -c $^

PoolDict.cxx:	$(LOCALINC)/Pool.hxx
				rootcint -f $@ -c $^

IsotopicVectorDict.cxx:		$(LOCALINC)/IsotopicVector.hxx
				rootcint -f $@ -c $^

ZAIDict.cxx:			$(LOCALINC)/ZAI.hxx
				rootcint -f $@ -c $^

FabricationPlantDict.cxx:	$(LOCALINC)/FabricationPlant.hxx
				rootcint -f $@ -c $^

EvolutionDataDict.cxx:	$(LOCALINC)/EvolutionData.hxx
				rootcint -f $@ -c $^




clean:
		@rm -vf $(OBJS) *~ core *Dict.cxx *Dict.h
oclean :
		@rm -vf $(OBJS) *~ core *Dict.cxx *Dict.h
.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx
