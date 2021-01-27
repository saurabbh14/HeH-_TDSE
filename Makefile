FC  		= gfortran 
FFLAGS		= -fopenmp -O3 -w 
LDFLAGS 	= -lfftw3 -lfftw3_omp -lm -llapack -lblas

default: 	dynamics

dynamics:   	main.o   adiabatic.o\
	propagation.o	setpot.o\
	prop2d.o	nuclear_wv.o
				
		${FC} *.o ${FFLAGS} ${LDFLAGS} -o dynamics
		
main.o:	 main.f90
	${FC} main.f90 ${FFLAGS} -c -o main.o 

setpot.o:	setpot.f90
	${FC} setpot.f90 ${FFLAGS} -c -o setpot.o 	
		
adiabatic.o:	adiabatic.f90
	${FC} adiabatic.f90 ${FFLAGS} -c -o adiabatic.o 
	
propagation.o:	propagation.f90
	${FC} propagation.f90 ${FFLAGS} -c -o propagation.o 
		
prop2d.o:	prop2d.f90
	${FC} prop2d.f90 ${FFLAGS} -c -o prop2d.o 
		
nuclear_wv.o:	nuclear_wv.f90
	${FC} nuclear_wv.f90 ${FFLAGS} -c -o nuclear_wv.o

clean:
	rm -f *.o 
	rm -f *.mod 
	rm -f dynamics 

