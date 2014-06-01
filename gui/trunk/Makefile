include ./Makefile.config

OMPFLAGS=
OMPLIB=

# Directory containing libraries
LIBDIR          = $(CLASS_lib)
# Directory containing includes for CLASS
LOCALINC	= $(CLASS_include)

######################################################
######### nothing to change from here ################
######################################################

#------------------------------------------------------------------------------
CXXFLAGS      += -I$(LOCALINC)
LD            = $(CXX)
LIBS          += -L$(LIBDIR) -lCLASSpkg -lgomp
LDFLAGS       = 
GLIBS         = $(ROOTGLIBS)

CXXFLAGS     += $(ROOTCFLAGS) 

#------------------------------------------------------------------------------
subO	= CLASSPlotElement.o \
	  CLASSWin.o \
	  CLASSRead.o
#------------------------------------------------------------------------------
subS	= CLASSPlotElement.cxx \
	  CLASSWin.cxx \
	  CLASSRead.cxx

#
mainO	=  CLASSGui.o
mainS	=  CLASSGui.cxx
main	=  CLASSGui


OBJS    = $(subO) $(mainO) 

all:    $(OBJS)
		
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -o $(main)
	mv  $(main) bin/
	@echo "$(main) done"

clean:
		@rm -f $(OBJS)

.SUFFIXES: .cxx

###
CLASSRead.o:	CLASSRead.hxx \
		CLASSPlotElement.hxx

CLASSWin.o:	CLASSRead.hxx \
		CLASSPlotElement.hxx \
		CLASSWin.hxx 

CLASSPlotElement.o:	CLASSPlotElement.hxx

.cxx.o:
	$(CXX) $(CXXFLAGS)  -c  $<

