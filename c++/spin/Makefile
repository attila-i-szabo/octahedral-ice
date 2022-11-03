# Makefile to compile the spinice simulator classes

CXX = icpc # g++
CXXFLAGS = -std=c++11 -pedantic -g -O3 -Wall

INC_PATH = -I ~/include -I ../common

LD = ld

# Lists all object files to be made for this project
all: spin.o octa.o spin_correlator.o ewald_rs.o ewald_generator.o spinice.o
# Clears the directory of object files
clean:
	rm *.o
# Recompiles all object files by first removing them
remake: clean all

# Nontrivial header dependences
spin.o: octa.hh
octa.o: spin.hh
ewald_generator.o: spin.hh
ewald_rs.o: ewald_generator.hh
spin_correlator.o: ewald_generator.hh
spinice.o: spin.hh octa.hh spin_correlator.hh ewald_rs.hh basic_stat.hh

# General .cc -> .o rule
%.o: %.cc %.hh 
	$(CXX) $(CXXFLAGS) $(INC_PATH) $< -c -o $@
