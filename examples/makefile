# directory paths
AMREXDIR = ../libamrex
VPATH = $(AMREXDIR)/lib

# compiler
CXX = g++
CPPFLAGS = -Wall -O3

# dependencies
HEADERS =
DEPLIBS = -lamrex
LIBS = -lgfortran

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -c $< -I$(AMREXDIR)/include -I../include

pulse: pulse.o $(DEPLIBS)
	$(CXX) $< -o $@ -L$(AMREXDIR)/lib $(DEPLIBS) $(LIBS)

all: pulse

.PHONY: clean

clean:
	rm -f *.o pulse
