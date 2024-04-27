OBJECTS = LinAl.o fgm.o
	
DOBJECTS = LinAl.o 

MODULES = LinAl.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas #-g3 -fbacktrace -fbounds-check

.PHONY: clean

output.txt: fgm.exe
	./fgm.exe > output.txt

fgm.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o fgm.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) fgm.exe *.dat
