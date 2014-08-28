# Directory containing libraries
LIBDIR          = ../../lib
# Base directory containing includes for libmctal and libvalerr
EXTERNINC	= ../external
# Directory containing includes for CLASS
LOCALINC	= ../include


EQM		= ../Model/Equivalence
IM		= ../Model/Irradiation
XSM		= ../Model/XS


ROOTCFLAGS    	= `root-config --cflags`
ROOTLIBS	= `root-config --libs`
ROOTGLIBS	= `root-config --glibs`

######### nothing to change from here #########
INCLUDES	= $(LOCALINC)/*.hxx 
LIBNAME 	= CLASSpkg
OBJS		= 	CLASSLogger.o \
			ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAIMass.o \
			CLASSObject.o CLASSObjectDict.o\
			CLASSFacility.o CLASSFacilityDict.o\
			FabricationPlant.o FabricationPlantDict.o \
			Reactor.o ReactorDict.o \
			CLASSBackEnd.o CLASSBackEndDict.o\
			Storage.o StorageDict.o\
			Pool.o PoolDict.o\
			DecayDataBank.o \
			DynamicalSystem.o\
			IrradiationModel.o \
			EquivalenceModel.o \
			XSModel.o \
			CLASSFuel.o\
			PhysicModels.o \
			EvolutionData.o EvolutionDataDict.o \
			CLASSFuelPlan.o\
			Scenario.o

OBJMODEL	=	$(EQM)/EQM_MLP_PWR_MOX.o $(EQM)/EQM_QUAD_PWR_MOX.o $(EQM)/EQM_LIN_PWR_MOX.o \
			$(XSM)/XSM_MLP.o $(XSM)/XSM_CLOSEST.o \
			$(IM)/IM_RK4.o $(IM)/IM_Matrix.o


ROOTOBJS	= 	CLASSLogger.o \
			ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAIMass.o \
			CLASSObject.o CLASSObjectDict.o\
			CLASSFacility.o CLASSFacilityDict.o\
			Reactor.o ReactorDict.o \
			FabricationPlant.o FabricationPlantDict.o \
			CLASSBackEnd.o CLASSBackEndDict.o\
			Storage.o StorageDict.o\
			Pool.o PoolDict.o\
			DecayDataBank.o \
			DynamicalSystem.o\
			CLASSFuel.o\
			EvolutionData.o EvolutionDataDict.o \
			CLASSFuelPlan.o\
			PhysicModels.o


CXX           = g++
CXXFLAGS      = -O2 -g -Wall -fopenmp -fPIC -I$(LOCALINC) $(ROOTCFLAGS)
LD            = g++
LDFLAGS       =  -g -Wall -fPIC $(ROOTLIBS) -lTMVA -shared -lgomp

all:		$(OBJS) $(OBJMODEL)
		$(LD) $(LDFLAGS) $(OBJS) $(OBJMODEL) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"
		$(LD) $(LDFLAGS) $(ROOTOBJS) -o $(LIBDIR)/lib$(LIBNAME)_root.so
		@echo "lib$(LIBNAME)_root.so done"

clean:
		@rm -vf $(OBJS) $(OBJMODEL) *~ core *Dict.cxx *Dict.h


install:	
		@ln -sf ../Model/* ../include/



CLASSObjectDict.cxx:		$(LOCALINC)/CLASSObject.hxx
				rootcint -f $@ -c $^

CLASSFacilityDict.cxx:		$(LOCALINC)/CLASSFacility.hxx
				rootcint -f $@ -c $^

CLASSBackEndDict.cxx:		$(LOCALINC)/CLASSBackEnd.hxx
				rootcint -f $@ -c $^

StorageDict.cxx:		$(LOCALINC)/Storage.hxx
				rootcint -f $@ -c $^

ReactorDict.cxx:		$(LOCALINC)/Reactor.hxx
				rootcint -f $@ -c $^
				
FabricationPlantDict.cxx:	$(LOCALINC)/FabricationPlant.hxx
				rootcint -f $@ -c $^

PoolDict.cxx:			$(LOCALINC)/Pool.hxx
				rootcint -f $@ -c $^

IsotopicVectorDict.cxx:		$(LOCALINC)/IsotopicVector.hxx
				rootcint -f $@ -c $^

ZAIDict.cxx:			$(LOCALINC)/ZAI.hxx
				rootcint -f $@ -c $^

EvolutionDataDict.cxx:		$(LOCALINC)/EvolutionData.hxx
				rootcint -f $@ -c $^


.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx -o $*.o
