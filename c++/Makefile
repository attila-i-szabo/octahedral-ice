# Makefile to compile and link with octahedral spin ice code

CXX = icpc # g++
CXXFLAGS = -std=c++11 -pedantic -g -O3 -Wall
LDFLAGS = -Wl,--gc-sections

LIBS_STD = -lstdc++ -lm
# LIBS_GSL = -lgsl -lgslcblas
LIBS_FFT = -lfftw3

SPIN_DIR = spin
SPIN_CC = $(wildcard $(SPIN_DIR)/*.cc)
SPIN_O = $(patsubst %.cc,%.o,$(SPIN_CC))

INC_PATH = -I $(SPIN_DIR) -I common

# Nontrivial header dependences
dipolar_correlators.o exchange_correlators.o otsuka_correlators.o: spin/spinice.hh spin/basic_stat.hh spin/ewald_rs.hh spin/ewald_generator.hh spin/spin_correlator.hh common/vec3.hh
all_samples.o: spin/spinice.hh spin/ewald_rs.hh spin/ewald_generator.hh common/vec3.hh

# Generic .cc -> .o rule: note we generally don't have .hh files
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INC_PATH) $< -c -o $@

# Make an object file in another directory
$(SPIN_DIR)/%.o:
	$(MAKE) -C $(dir $@) $(notdir $@)

# Nontrivial object dependences
dipolar_correlators.out exchange_correlators.out otsuka_correlators.out all_samples.out: $(SPIN_O)

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
# explicitly: 'library' object files aren't linked unless needed
%.out: %.o
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(LDFLAGS) $^ $(LIBS_STD) $(LIBS_FFT) \
	-o $@

# Don't allow direct compilation 
%.out: %.cc
# Don't allow removal of object files in chained implicit linking
.SECONDARY:
