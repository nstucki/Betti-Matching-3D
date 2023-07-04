# Update the following variables for your build setup.
CXX=g++
CPPFLAGS=-I /usr/local/include
LDFLAGS+=-L/usr/local/lib

# NO CONFIGURATION NEEDED BELOW THIS LINE
# =======================================

SRC=\
	minimize.cpp factor.cpp chunk.cpp reductions.cpp bireductions.cpp \
	lw.cpp complexes.cpp matrices.cpp ArrayColumn.cpp HeapColumn.cpp \
	time_measurement.cpp block_column_matrix.cpp Cone.cpp \
	relative_cohomology.cpp computation.cpp


FLAGS=-fopenmp #-fdiagnostics-color
FLAGS.release=-O3 -flto
FLAGS.debug=-g

CPPFLAGS+=$(FLAGS) -std=c++2a -Wall
CPPFLAGS.release=$(FLAGS.release) -DNDEBUG
CPPFLAGS.debug=$(FLAGS.debug) -D_LIBCPP_ENABLE_ASSERTIONS=1

LDFLAGS+=$(FLAGS)
LDFLAGS.release=$(FLAGS.release)
LDFLAGS.debug=$(FLAGS.debug)
LDLIBS=-lboost_program_options

DEP=dep
CPPFLAGS+=-MT $@ -MMD -MP -MF ${@F:%.o=${DEP}/%.d}


# default target
all: 2pac

# delete implicit rules
.SUFFIXES:

release/%: BUILD=release
debug/%: BUILD=debug

.SECONDEXPANSION:

# Linker targets
%/2pac: %/2pac.o $(patsubst %.cpp,$$*/%.o,$(SRC))
	$(CXX) $(LDFLAGS) $(LDFLAGS.$(BUILD))  -o $@ $^ $(LDLIBS)

# Compiler targets
%.o: $$(*F).cpp
	@mkdir -p $(@D) ${DEP}
	$(CXX) $(CPPFLAGS) $(CPPFLAGS.$(BUILD)) -o $@ -c $<

.PHONY: clean compile_commands.json
clean:
	rm -rf release debug ${DEP}

compile_commands.json: Makefile
	bear -- make --always make

DEPFILES=${SRC:%.cpp=${DEP}/%.d}
include $(wildcard $(DEPFILES))