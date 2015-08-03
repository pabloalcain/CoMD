# -*- Makefile -*-
SHELL=/bin/sh
CXX=g++
CFLAGS=-Wall -g -O0 -pg -ffast-math -Wextra -pedantic

# list of source files
SRC=box.cpp integrator.cpp particles.cpp potential.cpp sys.cpp cell.cpp dump.cpp comd.cpp 
OBJ=$(SRC:%.cpp=%.o)

default: comd.x

clean:
	rm -f *.o comd.x

# linker rule
comd.x: $(OBJ)
	$(CXX) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.cpp %.h
	$(CXX) -c $(CFLAGS) $<
