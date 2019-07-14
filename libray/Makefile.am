CFLAGS=-Wall -g
AROPT=

COMPILE_FLAG=@COMPILE_FLAG@
LINK_FLAG=@LINK_FLAG@

inc=-I../include/

OBJS=ray.o build.o solve.o ray_dump.o ray_stats.o ray_mat.o

all: $(OBJS)
	$(AR) rcs libgpvdm_ray.@SLIB_EXT@ $(OBJS)
	ln -s -f ./libray/libgpvdm_ray.@SLIB_EXT@ ../libgpvdm_ray.@SLIB_EXT@
#-cvq

.c.o : 
	$(CC) -c $(lib_type) $(CFLAGS) $(COMPILE_FLAG) $< -o $@  $(inc)

clean:
	rm *.o *.dll *.so *.a *.lib -f
