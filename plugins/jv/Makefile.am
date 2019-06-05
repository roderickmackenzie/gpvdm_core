OBJS=plugin.o jv.o

lib=-L../../
inc=-I../../include/

all: $(OBJS)
	$(CC) -shared $(OBJS) @DLL_LINKER@ $(lib) -o plugin.@DLIB_EXT@
	cp plugin.@DLIB_EXT@ ../jv.@DLIB_EXT@ -f

.c.o : 
	$(CC) $(inc) @DLL_COMPILE@ @COMPILE_FLAG@ -Wno-error=unused-but-set-variable -Wno-error=unused-variable -s -c $< -o $@

clean:
	rm *.o *.dll *.so -f

