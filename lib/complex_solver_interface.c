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

/** @file complex_solver_interface.c
	@brief Load the complex sparse matrix solver .so/.dll.  If this is UMFPACK the plugin will call UMFPACK, for other custom solvers the work will be done in the plugin.
*/


#include <stdio.h>
#include <stdlib.h>

	#include <dlfcn.h>

#include "util.h"
#include "inp.h"
#include "gpvdm_const.h"
#include "device.h"
#include "dump_ctrl.h"
#include "config.h"
#include "cal_path.h"
#include <lang.h>
#include <log.h>

static int unused __attribute__((unused));

void test_complex_solver_init(struct simulation *sim,char *solver_name)
{
char lib_path[PATH_MAX];

find_dll(sim, lib_path,solver_name);

char *error;

	sim->dll_complex_matrix_handle = dlopen(lib_path, RTLD_LAZY |RTLD_GLOBAL);

	if (!sim->dll_complex_matrix_handle)
	{
		ewe(sim, "%s\n", dlerror());
	}

	sim->dll_complex_matrix_solve = dlsym(sim->dll_complex_matrix_handle, "dll_complex_matrix_solve");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	sim->dll_complex_matrix_solver_free = dlsym(sim->dll_complex_matrix_handle, "dll_complex_matrix_solver_free");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	sim->dll_complex_matrix_init = dlsym(sim->dll_complex_matrix_handle, "dll_complex_matrix_init");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}


	(*sim->dll_complex_matrix_init)(sim);
}


void test_complex_solver_free(struct simulation *sim)
{
if (sim->dll_complex_matrix_handle!=NULL)
{
	(*sim->dll_complex_matrix_solver_free)(sim);

	printf_log(sim,"%s=%p\n",_("Freeing memory"),sim->dll_complex_matrix_handle);

	if (dlclose(sim->dll_complex_matrix_handle)!=0)
	{
		ewe(sim,"%s\n",_("Error closing dll"));
	}

	sim->dll_complex_matrix_handle=NULL;
	sim->dll_complex_matrix_solve=NULL;
	sim->dll_complex_matrix_solver_free=NULL;

}
}

