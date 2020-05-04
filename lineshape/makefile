OBJDIR = obj
EXEC = test
vpath %.o $(OBJDIR)/
vpath %.c $(SRCDIR)/

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
	CFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(LIBRARIES)
	CC = gcc
endif
ifeq ($(UNAME), Linux)
	LIBRARIES = -lm -lgsl -lgslcblas -lfftw3 -lmpi -L/usr/lib
	INCLUDES = -I/usr/include
	CFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -O3 -std=c99 $(DEBUGFLAG) $(LIBRARIES)
	CC = gcc
endif

SOURCES = ./src/parameters.c \
	  ./src/functions.c \
	  ./src/specDens.c

VPATH := $(sort $(dir $(SOURCES)))
OBJECTS = $(patsubst %.c, $(OBJDIR)/%.o, $(notdir $(SOURCES)))

# the "common" object files
$(OBJDIR)/%.o : %.c
	@echo "creating $@ ..."
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJECTS)
	@echo "building output ..."
	$(CC) -o $(EXEC) $(OBJECTS) $(LFLAGS)


.PHONY: clean cleaner
clean:
	\rm -f $(OBJECTS)

cleaner:
	\rm -f $(EXEC) $(OBJECTS)
