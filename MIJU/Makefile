#BIN = ../../bin
BIN = ./
program = $(BIN)/multiOPT.x 
include ./config.h
objs =   cales_discon.o sub.inverse.o calmat.o multipleSourcesOPT2D_discons.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.optrpt  
