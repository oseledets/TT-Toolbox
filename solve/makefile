include Makefile.mex

all: solve3d solve3d_2

solve3d: solve3d.c
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims solve3d.c -lmwlapack -lmwblas $(LIBS2) FC=$(FC) FOPTIMFLAGS="$(FOPTIMFLAGS)" CC=$(CC) COPTIMFLAGS="$(COPTIMFLAGS)" CFLAGS="$(CFLAGS)" LD=$(CC); exit;'

solve3d_2: solve3d_2.c
	$(MATLAB) -nodisplay -r 'mex -v -largeArrayDims solve3d_2.c ../fmex/tt-fort/mytt.a -lmwlapack -lmwblas $(LIBS2) FC=$(FC)  FOPTIMFLAGS="$(FOPTIMFLAGS)" CC=$(CC) COPTIMFLAGS="$(COPTIMFLAGS)" CFLAGS="$(CFLAGS)" LD="$(LDR) "; exit;'


# all:
# 	 icc -c  -I/opt/local/matlab2011b/extern/include -DMATLAB_MEX_FILE  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -O3   "solve3d.c"
# 	 icc -O3 -pthread -shared -Wl,--version-script,/opt/local/matlab2011b/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "solve3d.mexa64"  solve3d.o  -lmwlapack -lmwblas -Wl,-rpath-link,/opt/local/matlab2011b/bin/glnxa64 -L/opt/local/matlab2011b/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++ -lrt
