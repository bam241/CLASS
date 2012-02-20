OBJS          = fenetre.o litm2.o
PROG          = litm2

LIBDIR        = ../../../lib
INCDIR        = ..

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -c -O -g -Wall -fPIC $(ROOTCFLAGS) -I$(INCDIR)
LINK          = g++
LINKFLAGS     = -O -g -Wall -fPIC $(ROOTGLIBS)  -L$(LIBDIR)

all:	$(OBJS)
	$(LINK) $(LINKFLAGS) $(OBJS)  -lvalerr -lmctal -lmyroot -o $(PROG)

clean:
	rm -f $(OBJS) $(PROG) *~ core

.SUFFIXES: .cxx

.cxx.o: $*.hxx
	$(CXX) $(CXXFLAGS) -c $*.cxx
