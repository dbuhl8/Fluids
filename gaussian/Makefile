

gaussian: 
		cd fortran; gfortran gaussian_mod.f90 main.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o emain.exe && \
		./emain.exe; rm emain.exe; cd ..

inv: 
		cd fortran; gfortran gaussian_mod.f90 inv.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o einv.exe && \
		./einv.exe; rm einv.exe; cd ..

gen:
		cd fortran;gfortran gaussian_mod.f90 gen.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o egen.exe && \
		./egen.exe; rm egen.exe;cd ..

partial: 
		cd fortran; gfortran gaussian_mod.f90 partial.f90 -L/usr/local/lib -llapack -L/usr/lib -lblas -o epartial.exe && \
		./epartial.exe; rm epartial.exe; cd ..

clean: 
		cd ~/Documents/Fluids/gaussian/fortran; rm *.exe

