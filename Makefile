CC = g++
FLAGS = -g -c --std=c++17

SOURCEDIR = src_2
BUILDDIR = build_2

EXECUTABLE = BettiMatching
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
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