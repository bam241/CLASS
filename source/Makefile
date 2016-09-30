#################### CLASS main Makefile #########################
# include config file
include ../config/Makefile.config

# configuration :
	# INC_dir : directory containing includes for CLASS
	# SRC_dir : directory containing sources of majorities classes of CLASS
	# EQM_dir : directory containing sources of Equivalence Models
	# IM_dir  : directory containing sources of Irradiation Models
	# XQM_dir : directory containing sources of XS Models
	# EXT_dir : directory containing sources of external classes of CLASS (e.g. Graph, StringLines)
	# OBJ_dir : directory containing object files (*.o) after compilation of sources

INC_dir = include
SRC_dir = src

EQM_dir = Model/Equivalence
IM_dir  = Model/Irradiation
XSM_dir = Model/XS

EXT_dir = external

######################## nothing to change from here #########################
OBJ_dir = obj
LIBNAME = CLASSpkg

# files in $SRC_dir and $INC_dir
SRC_files = $(wildcard $(SRC_dir)/*.cxx)
INC_files = $(wildcard $(INC_dir)/*.hxx)

# files need a dictionary by ROOT
OBJ_Dict  = ZAIDict.o \
			IsotopicVectorDict.o \
			CLASSObjectDict.o \
			FabricationPlantDict.o \
			ReactorDict.o \
			CLASSBackEndDict.o \
			CLASSFacilityDict.o \
			SeparationPlantDict.o \
			StorageDict.o \
			PoolDict.o \
			EvolutionDataDict.o

# all object files of CLASS (build from source and dictionry list)
OBJ  = $(SRC_files:$(SRC_dir)/%.cxx=$(OBJ_dir)/%.o)
OBJ += $(OBJ_Dict:%.o=$(OBJ_dir)/%.o)

# source files in Model/Equivalence, Model/Irradiation and Model/XS
EQM_files = $(wildcard $(EQM_dir)/*.cxx)
IM_files  = $(wildcard $(IM_dir)/*.cxx)
XSM_files = $(wildcard $(XSM_dir)/*.cxx)

# all object model files
OBJ_model = $(EQM_files:$(EQM_dir)/%.cxx=$(OBJ_dir)/%.o) \
			$(IM_files:$(IM_dir)/%.cxx=$(OBJ_dir)/%.o) \
			$(XSM_files:$(XSM_dir)/%.cxx=$(OBJ_dir)/%.o)

# source files in external
EXT_files = $(wildcard $(EXT_dir)/*.cxx)

# all object external files
OBJ_ext   = $(EXT_files:$(EXT_dir)/%.cxx=$(OBJ_dir)/%.o)

# root object files, I don't understand this one
ROOTOBJS_ = CLASSLogger.o \
			ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			ZAIMass.o ZAIHeat.o ZAITox.o \
			CLASSNucleiFiliation.o \
			CLASSObject.o CLASSObjectDict.o \
			CLASSFacility.o CLASSFacilityDict.o \
			Reactor.o ReactorDict.o \
			FabricationPlant.o FabricationPlantDict.o \
			CLASSBackEnd.o CLASSBackEndDict.o \
			Storage.o StorageDict.o \
			Pool.o PoolDict.o \
			SeparationPlant.o SeparationPlantDict.o \
			DecayDataBank.o \
			IrradiationModel.o \
			EvolutionData.o EvolutionDataDict.o \
			ReactorScheduler.o\
			PhysicsModels.o

ROOTOBJS  = $(ROOTOBJS_:%.o=$(OBJ_dir)/%.o)

CXXFLAGS += -I$(INC_dir) $(ROOTCFLAGS)
LD        = $(CXX)
LIBS      = -L$(LIBDIR) $(ROOTLIBS) -lTMVA $(OMPLIB)
LDFLAGS  += -shared

all: $(OBJ) $(OBJ_model) $(OBJ_ext)
	@echo "\033[33m"lib$(LIBNAME).so : "\033[0m"
	$(LD) $(LDFLAGS) $(OBJ) $(OBJ_model) $(OBJ_ext) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME).so
	@echo -e "\033[1m    lib$(LIBNAME).so done\033[0m"
	@echo "\033[33m"lib$(LIBNAME)_root.so : "\033[0m"
	$(LD) $(LDFLAGS) $(ROOTOBJS) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME)_root.so
	@echo -e "\033[1m    lib$(LIBNAME)_root.so done\033[0m"

clean:
	@echo -e "\033[41;97;1mclean\033[0m"
	rm -f $(OBJ) $(OBJ_model) $(OBJ_ext) $(SRC_dir)/*~ $(SRC_dir)/*Dict.cxx $(SRC_dir)/*Dict.h $(EXT_dir)/*Dict.cxx $(EXT_dir)/*Dict.h
	rm -f $(INC_dir)/$(EQM_dir) $(INC_dir)/$(IM_dir) $(INC_dir)/$(XSM_dir) $(INC_dir)/$(EXT_dir)

external_link:
	@ln -sf $$(pwd)/external/ $(INC_dir)/

install:	
	@ln -sf $$(pwd)/Model/* $(INC_dir)/

# build all dictionnary
$(SRC_dir)/%Dict.cxx: $(INC_dir)/%.hxx
	@echo "\033[36m"$(notdir $@) : "\033[0m"
	cd $(SRC_dir) && rootcint -f $(notdir $@) -c ../$^

$(EXT_dir)/%Dict.cxx: $(EXT_dir)/%.hxx
	@echo "\033[36m"$(notdir $@) : "\033[0m"
	cd $(EXT_dir) && rootcint -f $(notdir $@) -c ../$^


.SUFFIXES: .cxx

# build standard sources
$(OBJ_dir)/%.o:		$(SRC_dir)/%.cxx $(INC_files)
	@echo "\033[95m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# build EQM sources
$(OBJ_dir)/%.o:		$(EQM_dir)/%.cxx $(EQM_dir)/%.hxx $(INC_files)
	@echo "\033[35m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# build IM sources
$(OBJ_dir)/%.o:		$(IM_dir)/%.cxx $(IM_dir)/%.hxx $(INC_files)
	@echo "\033[92m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# build XSM sources
$(OBJ_dir)/%.o:		$(XSM_dir)/%.cxx $(XSM_dir)/%.hxx $(INC_files)
	@echo "\033[34m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# build external sources
$(OBJ_dir)/%.o:		$(EXT_dir)/%.cxx $(EXT_dir)/%.hxx $(INC_files)
	@echo "\033[91m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

#$(OBJ_dir)/GraphDict.o:	$(EXT_dir)/GraphDict.cxx $(EXT_dir)/Graph.hxx $(INC_dir)
#	@echo "\033[92m"$(notdir $@) : "\033[0m"
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
