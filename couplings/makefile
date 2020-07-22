GF = gfortran
SOURCES = coupling_calc.f08
OBJECTS = $(patsubst %.f08, %.o, $(notdir $(SOURCES)))
CFLAGS = -std=f2008 -ffree-form -Wall -fcheck=bounds
LFLAGS = -llapack
EXEC = coupling_calc

DEBUG = 0
ifeq ($(DEBUG), 1)
	DEBUGFLAGS = -g
else
	DEBUGFLAGS = -O2
endif

%.o : %.f08
	$(GF) $(DEBUGFLAGS) $(CFLAGS) -c $<

$(EXEC): $(OBJECTS)
	$(GF) $(OBJECTS) $(LFLAGS) -o $(EXEC)

.PHONY: cleaner
cleaner:
	\rm -f $(EXEC)
