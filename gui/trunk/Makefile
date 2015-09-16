include ../config/Makefile.config

OMPFLAGS=
OMPLIB=

# Directory containing includes for CLASS
LOCALINC	= ../source/include

######################################################
######### nothing to change from here ################
######################################################

#------------------------------------------------------------------------------
CXXFLAGS      += -g -I$(LOCALINC)
LD            = $(CXX)
LIBS          += -L$(LIBDIR) -lCLASSpkg $(OMPLIB)
LDFLAGS       =
GLIBS         = $(ROOTGLIBS) -lTMVA

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
	@mkdir -p bin
	@mv  $(main) bin/
	@echo "$(main) done"
	@echo "CLASSGui is now available in $LIBDIR/../gui/bin "



clean:
		@rm -f $(OBJS)
		@rm -f ./bin/CLASSGui


.SUFFIXES: .cxx

CLASSRead.o:	CLASSRead.hxx \
		CLASSPlotElement.hxx

CLASSWin.o:	CLASSRead.hxx \
		CLASSPlotElement.hxx \
		CLASSWin.hxx 

CLASSPlotElement.o:	CLASSPlotElement.hxx

.cxx.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c  $<

