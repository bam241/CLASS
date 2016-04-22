#################### CLASS main Makefile #########################
# include config file

include ../../config/Makefile.config
# Directory containing includes for CLASS
LOCALINC = ../include

EQMINC	 = ../Model/Equivalence
IMINC	 = ../Model/Irradiation
XSMINC	 = ../Model/XS


######### nothing to change from here #########
INCLUDES = $(LOCALINC)/*.hxx 
LIBNAME  = CLASSpkg
OBJS	 = CLASSLogger.o \
			ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAIMass.o ZAIHeat.o ZAITox.o \
			CLASSNucleiFiliation.o \
			CLASSObject.o CLASSObjectDict.o\
			CLASSFacility.o CLASSFacilityDict.o\
			FabricationPlant.o FabricationPlantDict.o \
			Reactor.o ReactorDict.o \
			CLASSBackEnd.o CLASSBackEndDict.o\
			SeparationPlant.o SeparationPlantDict.o\
			Storage.o StorageDict.o\
			Pool.o PoolDict.o\
			DecayDataBank.o \
			DynamicalSystem.o\
			IrradiationModel.o \
			EquivalenceModel.o \
			XSModel.o \
			CLASSFuel.o\
			PhysicsModels.o \
			EvolutionData.o EvolutionDataDict.o \
			CLASSFuelPlan.o\
			Scenario.o

OBJMODEL = 	$(EQMINC)/EQM_PWR_MLP_MOX.o $(EQMINC)/EQM_PWR_MLP_MOX_Am.o \
			$(EQMINC)/EQM_PWR_POL_UO2.o\
			$(EQMINC)/EQM_FBR_BakerRoss_MOX.o $(EQMINC)/EQM_MLP_Kinf.o\
			$(EQMINC)/EQM_FBR_MLP_Keff.o \
			$(XSMINC)/XSM_MLP.o $(XSMINC)/XSM_CLOSEST.o \
			$(IMINC)/IM_RK4.o $(IMINC)/IM_Matrix.o


ROOTOBJS = 	CLASSLogger.o \
			ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAIMass.o ZAIHeat.o ZAITox.o \
			CLASSNucleiFiliation.o \
			CLASSObject.o CLASSObjectDict.o\
			CLASSFacility.o CLASSFacilityDict.o\
			Reactor.o ReactorDict.o \
			FabricationPlant.o FabricationPlantDict.o \
			CLASSBackEnd.o CLASSBackEndDict.o\
			Storage.o StorageDict.o\
			Pool.o PoolDict.o\
			SeparationPlant.o SeparationPlantDict.o\
			DecayDataBank.o \
			IrradiationModel.o \
			CLASSFuel.o\
			EvolutionData.o EvolutionDataDict.o \
			CLASSFuelPlan.o\
			PhysicsModels.o


CXXFLAGS      += -I$(LOCALINC) $(ROOTCFLAGS)
LD             = $(CXX)
LIBS           = -L$(LIBDIR) $(ROOTLIBS) -lTMVA $(OMPLIB)
LDFLAGS       += -shared 

all:		$(OBJS) $(OBJMODEL)
		$(LD) $(LDFLAGS) $(OBJS) $(OBJMODEL) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"
		$(LD) $(LDFLAGS) $(ROOTOBJS) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME)_root.so
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

SeparationPlantDict.cxx:	$(LOCALINC)/SeparationPlant.hxx
				rootcint -f $@ -c $^

IsotopicVectorDict.cxx:		$(LOCALINC)/IsotopicVector.hxx
				rootcint -f $@ -c $^

ZAIDict.cxx:			$(LOCALINC)/ZAI.hxx
				rootcint -f $@ -c $^

EvolutionDataDict.cxx:		$(LOCALINC)/EvolutionData.hxx
				rootcint -f $@ -c $^


.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $*.cxx -o $*.o
