.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

include Makefile.in


all: mytt.a fort_mvk4_mex fort_amen_solve_mex uchol solve3d solve3d_2 tt_read tt_write

mytt.a:
	cd tt-fort; make

fort_mvk4_mex: fort_mvk4_mex.F90 mytt.a
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims fort_mvk4_mex.F90 tt-fort/mytt.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

fort_amen_solve_mex: mytt.a fort_amen_solve_mex.F90
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims fort_amen_solve_mex.F90 tt-fort/mytt.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'


uchol_full.o: uchol_full.c
	$(CC) $(COPT) -c uchol_full.c -o uchol_full.o

uchol: uchol.c uchol_full.o
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims uchol.c uchol_full.o -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'


solve3d: solve3d.c
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims solve3d.c -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

solve3d_2: solve3d_2.c mytt.a
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims solve3d_2.c tt-fort/mytt.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

eig3d_primme: eig3d_primme.c tt-fort/primme/primme.a
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims eig3d_primme.c tt-fort/primme/primme.a -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

tt_read: tt_read.c
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims tt_read.c -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

tt_write: tt_write.c
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims tt_write.c -lmwlapack -lmwblas $(MLLIB2) FC=$(FC) FOPTIMFLAGS="$(FOPT)"  CC=$(CC) COPTIMFLAGS="$(COPT) $(COPTIMFLAGS)" CFLAGS="$(MCFLAGS)" LD="$(LDR)"; exit;'

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
		rm -f *.o
