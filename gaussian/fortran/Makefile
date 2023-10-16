

gaussian: 
		gfortran gaussian_mod.f90 main.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o emain.exe;\
		./emain.exe;

inv: 
		gfortran gaussian_mod.f90 inv.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o einv.exe;\
		./einv.exe

gen:
		gfortran gaussian_mod.f90 gen.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o egen.exe;\
		./egen.exe;

partial: 
		gfortran gaussian_mod.f90 partial.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o epartial.exe;\
		./epartial.exe

clean: 
		rm *.exe

