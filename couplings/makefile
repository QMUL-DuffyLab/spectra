GF = gfortran
SOURCE = coupling_calc.f08
CFLAGS = -std=f2008 -ffree-form -Wall -fcheck=bounds
LFLAGS = -llapack
EXEC = coupling_calc

DEBUG = 0
ifeq ($(DEBUG), 1)
	DEBUGFLAGS = -g
else
	DEBUGFLAGS = -O2
endif

$(EXEC):
	$(GF) $(DEBUGFLAGS) $(SOURCE) $(LFLAGS) $(CFLAGS) -o $(EXEC)

.PHONY: clean
clean:
	\rm -f $(EXEC)
