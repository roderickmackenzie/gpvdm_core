//
// General-purpose Photovoltaic Device Model gpvdm.com- a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
//
// Copyright (C) 2012-2017 Roderick C. I. MacKenzie info at gpvdm dot com
//
// https://www.gpvdm.com
//
//
// This program is free software; you can redistribute it and/or modify it
// under the terms and conditions of the GNU Lesser General Public License,
// version 2.1, as published by the Free Software Foundation.
//
// This program is distributed in the hope it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
//

/** @file newton_interface.c
	@brief Load and run the newton solve .dll/.so file.  They are hot swappable hence the interface.
*/



#include <stdio.h>
#include <stdlib.h>
#include "util.h"

	#include <dlfcn.h>

#include "inp.h"
#include "light_interface.h"
#include "gpvdm_const.h"
#include "device.h"
#include "dump_ctrl.h"
#include "config.h"
#include "cal_path.h"
#include "lang.h"
#include "log.h"
#include "newton_interface.h"

static int unused __attribute__((unused));


void newton_init(struct simulation *sim,char *solver_name)
{
//printf_log(sim,_("Solver initialization\n"));
char lib_path[1000];

find_dll(sim, lib_path,solver_name);


char *error;

	sim->dll_solver_handle = dlopen(lib_path, RTLD_LAZY);

	if (!sim->dll_solver_handle)
	{
		ewe(sim,"%s\n", dlerror());
	}

	sim->dll_solve_cur = dlsym(sim->dll_solver_handle, "dll_solve_cur");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	sim->dll_solver_realloc = dlsym(sim->dll_solver_handle, "dll_solver_realloc");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	sim->dll_solver_free_memory = dlsym(sim->dll_solver_handle, "dll_solver_free_memory");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}



}


void newton_set_min_ittr(struct device *in,int ittr)
{
in->newton_min_itt=ittr;
}

void solver_realloc(struct simulation *sim,struct device * in)
{
(*sim->dll_solver_realloc)(sim,in);
}

void solver_free_memory(struct simulation *sim,struct device * in)
{
(*sim->dll_solver_free_memory)(sim,in);
}

void newton_interface_free(struct simulation *sim)
{
dlclose(sim->dll_solver_handle);
}
