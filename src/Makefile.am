OBJS=  sim_find_n0.o sim_run.o pos.o initial.o config.o anal.o mesh.o run_electrical_dll.o sim.o

AROPT=

inc=-I../include/ -L../

all: main.c $(OBJS)
	@WIN_ICON@
	$(CC) main.c $(OBJS) -o ../gpvdm_core@EXE_EXT@ @WIN_RES_FILE@ $(inc) @COMPILE_FLAG@ -Wl,--whole-archive @LOCAL_LINK@ -Wl,--no-whole-archive @LDFLAGS@ @LIBS@ @OPEN_CL_LINK@
	ln -s -f gpvdm_core@EXE_EXT@ ../go@EXE_EXT@
	ln -s -f ./gpvdm_core/gpvdm_core@EXE_EXT@ ../../go@EXE_EXT@
	@echo "$(DESTDIR)$(bindir)" > bindir
	@echo "$(DESTDIR)$(datadir)" > datadir
	@echo "$(DESTDIR)$(datarootdir)" > datarootdir

	cp ../../sim.gpvdm ../base.gpvdm -f
	rm ../../base.gpvdm -f
	ln -s ./gpvdm_core/base.gpvdm ../../base.gpvdm
	

.c.o : 
	$(CC) -c $(DEFINE_FLAGS) $(inc) $(CPPFLAGS) $(warn) @COMPILE_FLAG@ $< -o $@

clean:
	rm *.o *.a *.dll *.exe ../../base.gpvdm -f
