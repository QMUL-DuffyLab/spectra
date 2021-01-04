OBJDIR = obj
EXEC = ver_test
vpath %.o $(OBJDIR)/
vpath %.cpp $(SRCDIR)/

UNAME := $(shell uname)
ifneq ($(shell which git 2> /dev/null),)
	BRANCHNAME := $(shell git name-rev --name-only HEAD)
endif
ifneq ($(BRANCHNAME),)
	EXEC_NAME = $(EXEC)_$(UNAME)_$(BRANCHNAME)
else
	EXEC_NAME = $(EXEC)_$(UNAME)
endif

# Debug?
DEBUG = 1
ifeq ($(DEBUG), 1)
	DEBUGFLAG = -g
else
	DEBUGFLAG = -O2
endif

ifeq ($(UNAME), Darwin)
	LIBRARIES = -lm -lgsl -lfftw3 -L/opt/local/lib -L/opt/local/lib/mpich-devel-mp -L/opt/local/lib/openmpi-mp
	INCLUDES = -I/opt/local/include
	CFLAGS = -Wall -Wextra -std=gnu99 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -std=gnu99 $(DEBUGFLAG) $(LIBRARIES)
	CC = g++
endif
ifeq ($(UNAME), Linux)
	LIBRARIES = -lm -lgsl -lgslcblas -lfftw3 -llapacke -lmpi -L/usr/lib
	INCLUDES = -I/usr/include
	CFLAGS = -Wall -Wextra -std=c++11 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -std=c++11 $(DEBUGFLAG) $(LIBRARIES)
	CC = g++
endif

SOURCES = ./src/LSODA.cpp \
	  ./src/ver.cpp

VPATH := $(sort $(dir $(SOURCES)))
OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(notdir $(SOURCES)))

# the "common" object files
$(OBJDIR)/%.o : %.cpp
	@echo "creating $@ ..."
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJECTS)
	@echo "building output ..."
	$(CC) -o $(EXEC) $(OBJECTS) $(LFLAGS)


.PHONY: clean cleaner
clean:
	\rm -f $(OBJECTS)

cleaner:
	\rm -f $(EXEC) $(OBJECTS) $(DEP)
