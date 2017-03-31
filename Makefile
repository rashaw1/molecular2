# Project name
PROJECT = molecular.x
# Compiler
CXX = clang++

# Run Options
COMMANDLINE_OPTIONS = 

# Compiler options
DEBUG = -g -Wall -O0 -Wno-unused-variable -Wno-sign-compare -std=c++11 -D_GLIBCXX_DEBUG -omp -DLIBINT2_MAX_DERIV_ORDER=2 
OPTIM = -O3 -Wall -Wno-unused-variable -Wno-sign-compare -std=c++11 -omp -DLIBINT2_MAX_DERIV_ORDER=2
COMPILE_OPTIONS = $(OPTIM)

# Header include directories
HEADERS = -I./inc -I/usr/local/Cellar/eigen/3.3.1/include/eigen3 -I/usr/local/libint/2.3.0-beta.3/include -I/usr/local/libint/2.3.0-beta.3/include/libint2

# Libraries for linking
LIBS = -L/usr/local/libint/2.3.0-beta.3/lib -lint2 -std=c++11

# Dependency options
DEPENDENCY_OPTIONS = -MM

#-- Do not edit below this line --

# Subdirs to search for additional source files
SUBDIRS := $(shell ls -F | grep "\/")
DIRS := ./ $(SUBDIRS)
SOURCE_FILES := $(foreach d, $(DIRS), $(wildcard $(d)*.cpp) )

# Create an object for every cpp file
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES := $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Create .d files
%.d: %.cpp
	$(CXX) $(DEPENDENCY_OPTIONS) $< -MT "$*.o $*.d" -MF $*.d

# Make $(PROJECT) the default target
all: $(DEPENDENCIES) $(PROJECT)

$(PROJECT): $(OBJECTS)
	$(CXX) -o $(PROJECT) $(OBJECTS) $(LIBS)

# Include dependencies

# Compile every cpp file to an object
%.o: %.cpp
	$(CXX) -c $(COMPILE_OPTIONS) -o $@ $< $(HEADERS)

# Build and run project
run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

# Clean and debug
.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(PROJECT) $(DEPENDENCIES)

clean-all: clean depclean
