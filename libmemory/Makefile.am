inc=-I../include/
OBJS=memory.o memory_basic.o zxy_long_double.o matrix.o zxy_long_double_complex.o zy_long_double.o zxy_int.o light_l_long_double.o light_y_long_double.o light_zxyl_long_double.o light_zxyl_long_double_complex.o light_zxy_long_double.o zx_long_double.o light_zxy_p_object.o heat_zxy_long_double.o heat_zxy_p_object.o zx_epitaxy_int.o

COMPILE_FLAG=@COMPILE_FLAG@
LINK_FLAG=@LINK_FLAG@

all: $(OBJS)
	$(AR) cr libgpvdm_memory.@SLIB_EXT@ $(OBJS)
	ln -s -f ./libmemory/libgpvdm_memory.@SLIB_EXT@ ../libgpvdm_memory.@SLIB_EXT@
#-cvq

.c.o:
	$(CC) -c $(CPPFLAGS) $(COMPILE_FLAG) $(inc) $< -o $@  

clean:
	rm *.o *.dll *.so *.a *.lib -f
