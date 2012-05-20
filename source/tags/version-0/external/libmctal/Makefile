LIBDIR        = ../../../lib
INCDIR        = ..

# Nothing to change under this line


LIBNAME       = mctal
OBJS          = TMComment.o TMTFCData.o TMTFC.o TMTally.o TMctal.o TMKcode.o
INCLUDES      = *.hxx

CXX           = g++
CXXFLAGS      = -O -g -Wall -fPIC -I $(INCDIR)
LD            = g++
LDFLAGS       = -g -Wall -L $(LIBDIR) -lvalerr -shared

all:		$(OBJS)
		$(LD) $(LDFLAGS) $(OBJS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"

clean:
		@rm -f $(OBJS) *~ core

.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx
