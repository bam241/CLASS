ROOTCFLAGS:=$(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTGLIBS:=$(shell ${ROOTSYS}/bin/root-config --glibs)
ROOTLIBS:=$(shell ${ROOTSYS}/bin/root-config --libs)

CXX= g++ -std=c++11
CXXFLAGS= -O2 -g -fPIC -finline-functions

LIBDIR=lib/

SRC_dir = src/
INC_dir = inc/

SRC_files = $(wildcard $(SRC_dir)/*.cxx)
INC_files = $(wildcard $(INC_dir)/*.hxx)

OBJ_dir = obj/
LIBNAME = Runpkg

# files need a dictionary by ROOT
OBJ_Dict  = ZAIDict.o \
			IsotopicVectorDict.o \
			ReactorDict.o \
			FacilityDict.o \
			RunDict.o

# all object files (build from source and dictionry list)
OBJ  = $(SRC_files:$(SRC_dir)/%.cxx=$(OBJ_dir)/%.o)
OBJ += $(OBJ_Dict:%.o=$(OBJ_dir)/%.o)

# root object files, I don't understand this one
ROOTOBJS_ = ZAI.o ZAIDict.o \
			IsotopicVector.o IsotopicVectorDict.o \
			Reactor.o ReactorDict.o \
			Facility.o FacilityDict.o \
			Run.o RunDict.o

ROOTOBJS  = $(ROOTOBJS_:%.o=$(OBJ_dir)/%.o)

CXXFLAGS += -I$(INC_dir) $(ROOTCFLAGS)
LD        = $(CXX)
LIBS      = -L$(LIBDIR) $(ROOTLIBS) -lTMVA $(OMPLIB)
LDFLAGS  += -shared

all: $(OBJ) $(OBJ_model) $(OBJ_ext)
	@echo "\033[33m"lib$(LIBNAME).so : "\033[0m"
	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME).so
	@echo "\033[1m    lib$(LIBNAME).so done\033[0m"
	@echo "\033[33m"lib$(LIBNAME)_root.so : "\033[0m"
	$(LD) $(LDFLAGS) $(ROOTOBJS) $(LIBS) -o $(LIBDIR)/lib$(LIBNAME)_root.so
	@echo "\033[1m    lib$(LIBNAME)_root.so done\033[0m"

clean:
	@echo "\033[41;97;1mclean\033[0m"
	rm -f $(OBJ) $(OBJ_model) $(OBJ_ext) $(SRC_dir)/*~ $(SRC_dir)/*Dict.cxx $(SRC_dir)/*Dict.h

.SUFFIXES: .cxx

# build standard sources
$(OBJ_dir)/%.o:		$(SRC_dir)/%.cxx $(INC_files)
	@echo "\033[95m"$(notdir $@) : "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@
