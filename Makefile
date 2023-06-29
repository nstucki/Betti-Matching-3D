CC = g++
FLAGS = --std==c++17

SOURCEDIR = /src/code
BUILDDIR = build

EXECUTABLE = BettiMatching
SOURCES = $(wildcard src/code/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))

all: dir $(BUILDDIR)/$(EXECUTABLE)

dir:
    mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXECUTABLE): $(OBJECTS)
    $(CC) $^ -o $@

$(OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.cpp
    $(CC) $(FLAGS) $< -o $@

clean:
    rm -f $(BUILDDIR)/*o $(BUILDDIR)/$(EXECUTABLE)