OBJECTS = LinAl.o shearSolve.o
	
DOBJECTS = LinAl.o 

MODULES = LinAl.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas #-g3 -fbacktrace -fbounds-check

.PHONY: clean

output.txt: shearSolve.exe
	./shearSolve.exe > output.txt

shearSolve.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o shearSolve.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) shearSolve.exe *.dat
