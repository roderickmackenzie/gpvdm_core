COMPILE_FLAG=@COMPILE_FLAG@
LINK_FLAG=@LINK_FLAG@

inc=-I../include/
OBJS=pl.o probe.o emission.o

all: $(OBJS)
	$(AR) rcs libgpvdm_emission.@SLIB_EXT@ $(OBJS)
	ln -s -f ./libemission/libgpvdm_emission.@SLIB_EXT@ ../libgpvdm_emission.@SLIB_EXT@
#-cvq

.c.o:
	$(CC) -c $(lib_type) $(CPPFLAGS) $(COMPILE_FLAG) $< -o $@  $(inc)

clean:
	rm *.o *.dll *.so *.a *.lib -f
