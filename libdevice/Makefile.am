inc=-I../include/
OBJS=device.o build_scene.o scene_dump.o add_shape.o epitaxy.o epitaxy_shapes.o device_malloc.o device_free.o device_heat.o epitaxy_optical.o epitaxy_heat.o device_world_stats.o

COMPILE_FLAG=@COMPILE_FLAG@
LINK_FLAG=@LINK_FLAG@

all: $(OBJS)
	$(AR) cr libgpvdm_device.@SLIB_EXT@ $(OBJS)
	ln -s -f ./libdevice/libgpvdm_device.@SLIB_EXT@ ../libgpvdm_device.@SLIB_EXT@
#-cvq

.c.o:
	$(CC) -c $(CPPFLAGS) $(COMPILE_FLAG) $(inc) $< -o $@  

clean:
	rm *.o *.dll *.so *.a *.lib -f
