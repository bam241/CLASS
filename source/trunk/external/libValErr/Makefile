LIBDIR        = ../../../lib
INCDIR        = ..

# Nothing to change under this line


LIBNAME       = valerr
OBJS          = ValErr.o
INCLUDES      = *.hxx

CXX           = g++
CXXFLAGS      = -O -g -Wall -fPIC -I $(INCDIR)
LD            = g++
LDFLAGS       = -g  -Wall -shared

all:		$(OBJS)
		$(LD) $(LDFLAGS) $(OBJS) -o $(LIBDIR)/lib$(LIBNAME).so
		@echo "lib$(LIBNAME).so done"

clean:
		@rm -f $(OBJS) *~ core

.SUFFIXES: .cxx

%.o: %.cxx $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $*.cxx
