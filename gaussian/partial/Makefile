
#Directory Stucture for Makefile
#Plot Directory
PFOLD = ~/Documents/Fluids/gaussian/plotData
#Log File Directory
LFOLD = ~/Documents/Fluids/gaussian/logs


FC = gfortran
FFLAGS = -llapack -lblas
#What do LFLAGS do??????
#LFLAGS = -L/usr/local/lib -L/usr/lib

#Module File
OMOD = gaussian_mod.o
FMOD = ../gaussian_mod.f90

.PHONY : clean

main: main.dat
		gnuplot '$(PFOLD)/main.p'

gen: gen.dat
		gnuplot '$(PFOLD)/gen.p'

partial: partial.dat
		gnuplot '$(PFOLD)/partial.p'

%.dat : %.e
		./$< >> $(LFOLD)/$(subst .dat,.txt,$@)

%.e : %.o
		$(FC) -o $@ $(OMOD) $< $(FFLAGS)

%.o : %.f90
		$(FC) -c $(FMOD) $< $(FFLAGS)

clean: 
		rm -rf *.e *.o

