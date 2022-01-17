# Warrick's attempt at a new Makefile.
# Drawn from
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

#FC = g77
FC = gfortran

#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15 -fbounds-check
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15 -m32
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O2
#FFLAGS = -extend_source -fast # -fpe3
# FFLAGS = -extend_source -g
FFLAGS = -ffixed-line-length-none -fno-automatic

#ifort options
#FFLAGS = -e -fast
#FFLAGS = -extend_source -ftz -static -C -fast

ODIR=obj
SDIR=src

_OBJ = main.o pressi.o statef.o fdirac.o consts.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

VPATH = $(SDIR)

$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)

FD: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ FD
