#BIN = ../../bin
BIN = ../bin/
program = $(BIN)/multiOPT.x 
program2 = $(BIN)/modelgen.x
include ../etc/config.h
objs =   modules.o others.o vectors.o disconConfig.o seismicSources.o videos.o   cales.o cales_discon.o cales_free.o sub.inverse.o calmat.o OPT2D.o
objs2 = modelgenerator.o

all : $(program) $(program2)

.SUFFIXES: .o .f90

$(program2): $(objs2)
	$(FC)  -o $@ $(objs2) $(FFLAGS)

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(program2) $(objs) $(objs2) *.optrpt  

