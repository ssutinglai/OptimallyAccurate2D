F90= gfortran
D1=	analitic.o fft.o cmul.o samax.o
D2=	gauss1.o
D3=	 samax.o

NORM:	$(D3)
	$(F90)   $(D3) -o norm

GAUSS:	$(D2)
	$(F90)   $(D2) -o gauss

ANA:	$(D1)
	$(F90)   $(D1) -o analitic

.SUFFIXES: .o .f .f90

.f.o:
	$(F90)  $*.f
.f90.o:
	$(F90)  $*.f90
