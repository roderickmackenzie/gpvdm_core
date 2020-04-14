OBJS=solver.o lib.o
CFLAGS=-Wall
inc=-I../../include/ -I/home/rod/Downloads/superlu_5.2.1/SuperLU_5.2.1/SRC/ #last line is for debian compatibility
lib=-L../../

all: $(OBJS)
	$(CC) -shared $(OBJS) @DLL_LINKER@ $(lib) @DLL_DEBUG@ $(LDFLAGS) -lsuperlu -o plugin.@DLIB_EXT@
	mv plugin.@DLIB_EXT@ ../superlu.@DLIB_EXT@

.c.o : 
	$(CC) -Wall -g $(inc) $(CPPFLAGS) @DLL_COMPILE@ @COMPILE_FLAG@ -c $< -o $@

clean:
	rm *.o *.dll *.so -f
