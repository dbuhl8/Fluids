
#Directory Stucture for Makefile
#Plot Directory
PFOLD = ~/Documents/Fluids/gaussian/plotData
TFOLD = ~/Documents/Fluids/gaussian/plots
#Log File Directory
LFOLD = ~/Documents/Fluids/gaussian/logs


FC = gfortran
FFLAGS = -llapack -lblas
#What do LFLAGS do??????
#LFLAGS = -L/usr/local/lib -L/usr/lib

#Module File
OMOD = gaussian_mod.o
FMOD = ../src/gaussian_mod.f90

.PHONY : clean

partial: partial.dat
		gnuplot $(PFOLD)/partial.p
		mv $(TFOLD)/partial.png $(TFOLD)/trials/$(shell tail -n 1 $(LFOLD)/$@.txt)
		tail -n 22 $(LFOLD)/$@.txt
		head -82000 $(PFOLD)/pgen.dat | tail +80000 > $(PFOLD)/pgen4045.dat
		#gwenview $(TFOLD)/trials/$(shell tail -n 1 $(LFOLD)/$@.txt)

%.dat : %.e
		./$< >> $(LFOLD)/$(subst .dat,.txt,$@)

%.e : %.o
		$(FC) -o $@ $(OMOD) $< $(FFLAGS)

%.o : %.f90
		$(FC) -c $(FMOD) $< $(FFLAGS)

clean: 
		rm -rf *.e *.o


