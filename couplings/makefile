GF = gfortran
SOURCES = ./src/aux.f08 \
	  ./src/coupling_calc.f08
OBJECTS = $(patsubst %.f08, %.o, $(notdir $(SOURCES)))
MOD_DIR = mod
CFLAGS = -J$(MOD_DIR) -std=f2008 -ffree-form -Wall -fcheck=bounds
LFLAGS = -llapack
EXEC = coupling_calc

VPATH := $(sort $(dir $(SOURCES)))
$(shell mkdir -p $(MOD_DIR))

DEBUG = 0
ifeq ($(DEBUG), 1)
	DEBUGFLAGS = -g -fbounds-check -ffpe-trap=denormal
else
	DEBUGFLAGS = -O2
endif

%.o : %.f08
	$(GF) $(DEBUGFLAGS) $(CFLAGS) -c $<

$(EXEC): $(OBJECTS)
	$(GF) $(OBJECTS) $(LFLAGS) -o $(EXEC)

.PHONY: clean cleaner
clean:
	\rm -f $(EXEC)
cleaner:
	\rm -f $(EXEC) $(MOD_DIR)
