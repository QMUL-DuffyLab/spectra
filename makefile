OBJDIR = obj
SRCDIR = src
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
	LIBRARIES = -lm -lgsl -L/opt/local/lib -L/opt/local/lib/mpich-devel-mp -L/opt/local/lib/openmpi-mp
	INCLUDES = -I/opt/local/include
	CFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(LIBRARIES)
	CC = gcc
endif
ifeq ($(UNAME), Linux)
	LIBRARIES = -lm -lgsl -lgslcblas -lmpi -L/usr/local/lib
	INCLUDES = -I/usr/local/include
	CFLAGS = -Wall -Wextra -std=c99 $(DEBUGFLAG) $(INCLUDES)
	LFLAGS = -Wall -Wextra -O3 -std=c99 $(DEBUGFLAG) $(LIBRARIES)
	CC = mpicc
endif

SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJTMP := $(notdir $(SOURCES))
OBJTMP := $(OBJTMP:.c=.o)
OBJ := $(addprefix $(OBJDIR)/, $(OBJTMP))
DEP = $(OBJ:.o=.d)

print-%  : ; @echo $* = $($*)

# pull in existing dependency info
-include $(DEP)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LFLAGS) $(OBJ) -o $(EXEC)

# Automatically build a list of dependencies for each object file
# http://scottmcpeak.com/autodepend/autodepend.html
# list dependencies and send to a file with extension .d
# moves the gcc output to a temp file so we can change formatting
# substitutes whatever's before the colon for its stem + .o;
# which allows make to be run from a different directory than src/obj files
# strips everything before the colon, then gets rid of continuation backslashes
# Finally strips leading spaces left after the backslashes are cut out
# and deletes the temporary dependency files
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)/
	$(CC) -c $(CFLAGS) $(SRCDIR)/$*.c -o $(OBJDIR)/$*.o
	$(CC) -MM $(CFLAGS) $(SRCDIR)/$*.c > $(OBJDIR)/$*.d
	@mv -f $(OBJDIR)/$*.d $(OBJDIR)/$*.d.tmp
	@sed -e 's|.*:|$(OBJDIR)/$*.o:|' < $(OBJDIR)/$*.d.tmp > $(OBJDIR)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $(OBJDIR)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(OBJDIR)/$*.d
	@rm -f $(OBJDIR)/$*.d.tmp

.PHONY: clean cleaner
clean:
	\rm -f $(OBJ)

cleaner:
	\rm -f $(EXEC) $(OBJ) $(DEP)
