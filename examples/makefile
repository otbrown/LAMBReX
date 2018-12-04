# directory paths
AMREXDIR = ../libamrex
LAMBREXDIR = ../liblambrex
VPATH = $(AMREXDIR)/lib:$(LAMBREXDIR)

# compiler
CXX = g++
CPPFLAGS = -Wall -O3

# dependencies
HEADERS =
DEPLIBS = -lamrex -llambrex
LIBS = -lgfortran

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -c $< -I../include

pulse: pulse.o $(DEPLIBS)
	$(CXX) $< -o $@ -L$(AMREXDIR)/lib -L$(LAMBREXDIR) $(DEPLIBS) $(LIBS)

all: pulse

.PHONY: clean

clean:
	rm -f *.o pulse
