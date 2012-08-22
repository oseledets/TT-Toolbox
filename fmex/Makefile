.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

include Makefile.in


all: mytt.a fort_mvk4_mex fort_amr_solve_mex

mytt.a:
	cd tt-fort; make

fort_mvk4_mex: fort_mvk4_mex.F90 mytt.a
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims fort_mvk4_mex.F90 tt-fort/mytt.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

fort_amr_solve_mex: mytt.a fort_amr_solve_mex.F90
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims fort_amr_solve_mex.F90 tt-fort/mytt.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'


#
# .f.o:
# 		$(FC) $(FOPT) -c $<
# .f90.o:
# 		$(FC) $(FOPT) -c $<
# .F90.o:
# 		$(FC) $(FOPT) -c $<
# .c.o:
# 		$(CC) $(COPT) -c $<


clean:
		rm -f fort_mvk4_mex.mexa64 fort_amr_solve_mex.mexa64
