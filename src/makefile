# directory paths
AMREXDIR = ../libamrex
LAMBREXDIR = ../liblambrex
INCLUDEDIR = ../include
VPATH = $(AMREXDIR)/lib:$(LAMBREXDIR):$(INCLUDEDIR)

# compiler
CXX = g++
CPPFLAGS = -Wall -O3

# dependencies
HEADERS = lambrex.h
DEPLIBS = -lamrex
LIBS = -lgfortran
OBJS = lambrex.o Lamb.o

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -c $< -I. -I$(AMREXDIR)/include -I$(INCLUDEDIR)

liblambrex.a: $(OBJS)
	ar rcs liblambrex.a $(OBJS)

.PHONY: install uninstall clean distclean

install: liblambrex.a
	cp liblambrex.a $(LAMBREXDIR)/

uninstall:
	rm -f $(LAMBREXDIR)/liblambrex.a

clean:
	rm -f *.o *.a

distclean: clean uninstall
