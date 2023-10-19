

FFOLD = ~/Documents/Fluids/gaussian/fortran
PFOLD = ~/Documents/Fluids/gaussian/plotData

FC = gfortran
FFLAGS = -llpack -lblas
LFLAGS = -L/usr/local/lib -L/usr/lib


.PHONY : clean

gaussian: 
		cd $(FFOLD); $(FC) gaussian_mod.f90 main.f90 $(LFLAGS) $(FFLAGS) -o emain.exe && \
		./emain.exe; rm emain.exe; cd ..
		@gnuplot 'main.p'

inv: 
		cd fortran; gfortran gaussian_mod.f90 inv.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o einv.exe && \
		./einv.exe; rm einv.exe; cd ..

gen:
		cd $(FFOLD);$(FC) gaussian_mod.f90 gen.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o egen.exe && \
		./egen.exe; rm egen.exe;cd ..
		@gnuplot 'plotData/gen.p'


partial: 
		cd fortran; gfortran gaussian_mod.f90 partial.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o epartial.exe && \
		./epartial.exe; rm epartial.exe; cd ..
		@gnuplot 'partial.p'

clean: 
		rm -r *.exe *.o

