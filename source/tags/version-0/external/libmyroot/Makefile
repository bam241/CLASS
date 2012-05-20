LIBDIR        = ../../../lib
INCDIR        = ..

# Nothing to change under this line

LIBNAME       = myroot
OBJS          = THistoCalc.o THisto2Calc.o  NextName.o
INCLUDES      = *.hxx

CXX           = g++
CXXFLAGS      = -O -g -Wall `root-config --cflags` -fPIC
LD            = g++
LDFLAGS       = -g -Wall `root-config --glibs` -L$(LIBDIR)/
SOFLAGS       = -shared

# $(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -llapack -lblas -lvalerr -o $(LIBDIR)/lib$(LIBNAME).so

# $(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -lblas -lvalerr -o $(LIBDIR)/lib$(LIBNAME).so
all:		$(OBJS)
		 $(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS)  -lvalerr -o $(LIBDIR)/lib$(LIBNAME).so
		
		@echo "lib$(LIBNAME).so done"

clean:
		@rm -f $(OBJS)  lib$(LIBNAME).so *~ core

.SUFFIXES: .cxx

.cxx.o: $*.hxx
	$(CXX) $(CXXFLAGS) -I $(INCDIR) -c $*.cxx
