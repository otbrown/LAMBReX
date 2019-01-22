# directory paths
AMREXDIR = ../libamrex
LAMBREXDIR = ../liblambrex
INCLUDEDIR = ../include
VPATH = $(INCLUDEDIR):$(AMREXDIR)/lib:$(LAMBREXDIR)

# compiler
CXX = g++
CPPFLAGS = -Wall -O3

# dependencies
HEADERS = lambrex.h
DEPLIBS = -llambrex -lamrex
LIBS = -lgfortran

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -c $< -I$(INCLUDEDIR) -I$(AMREXDIR)/include

pulse: pulse.o $(DEPLIBS)
	$(CXX) $< -o $@ -L$(LAMBREXDIR) -L$(AMREXDIR)/lib $(DEPLIBS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o pulse
