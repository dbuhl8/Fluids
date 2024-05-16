
COMPILER = mpif90

OBJECTS = shearSolve.o
#DOBJECTS = LinAl.o 

MODULES = #LinAl.mod

FFLAGS = -lblas -llapack -g #-g3 -fbacktrace -fbounds-check

.PHONY: clean

shearSolve: $(MODULES) $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(FFLAGS) -o shearSolve

%.o: %.f90
	$(COMPILER) -c $(FFLAGS) $<

%.mod: %.f90
	$(COMPILER) -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) shearSolve
