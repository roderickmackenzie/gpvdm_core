OBJS=plugin.o newton.o
lib=-L../../
inc=-I../../include/

all: $(OBJS)
	$(CC) -shared $(OBJS) @DLL_LINKER@ $(lib) @DLL_DEBUG@ -o plugin.@DLIB_EXT@
	mv plugin.@DLIB_EXT@ ../newton_simple.@DLIB_EXT@

.c.o : 
	$(CC) -Wall -g $(inc) $(CPPFLAGS) -Wno-unused-variable -Wno-unused-but-set-variable @DLL_COMPILE@ @COMPILE_FLAG@ -c $< -o $@

clean:
	rm *.o *.dll *.so -f
