inc=-I../include/ -I/usr/include/suitesparse/ #last line is for debian compatibility
OBJS=advmath.o rand.o vec.o triangle_io.o

COMPILE_FLAG=@COMPILE_FLAG@
LINK_FLAG=@LINK_FLAG@

all: $(OBJS)
	$(AR) cr libgpvdm_basicmath.@SLIB_EXT@ $(OBJS)
	ln -s -f ./libbasicmath/libgpvdm_basicmath.@SLIB_EXT@ ../libgpvdm_basicmath.@SLIB_EXT@
#-cvq

.c.o:
	$(CC) -c $(CPPFLAGS) $(COMPILE_FLAG) $(inc) $< -o $@  

clean:
	rm *.o *.dll *.so *.a *.lib -f
