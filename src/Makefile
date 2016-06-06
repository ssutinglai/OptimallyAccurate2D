#BIN = ../../bin
BIN = ../bin/
program = $(BIN)/multiOPT.x 
include ../etc/config.h
objs =   cales_discon.o cales_free.o sub.inverse.o calmat.o OPT2D.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.optrpt  
