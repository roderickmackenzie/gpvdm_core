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

/** @file light_interface.c
	@brief Acts as an interface between the model and the light dlls for each optical model.
*/

	#include <dlfcn.h>

#include "util.h"
#include "inp.h"
#include "light_interface.h"
#include "gpvdm_const.h"
#include "device.h"
#include "dump_ctrl.h"
#include "config.h"
#include "cal_path.h"
#include "lang.h"
#include "log.h"
#include "sim.h"
#include "memory.h"
#include <light_fun.h>
#include <matrix.h>

static int unused __attribute__((unused));

void light_load_dlls(struct simulation *sim,struct light *li)
{
	char lib_path[PATH_MAX];
	char lib_name[100];

	printf_log(sim,"%s\n",_("Initializing optical model"));

	sprintf(lib_name,"light_%s",li->mode);
	find_dll(sim, lib_path,lib_name);


	char *error;

	li->lib_handle = dlopen(lib_path, RTLD_LAZY);

	if (!li->lib_handle)
	{
		ewe(sim, "%s\n", dlerror());
	}

	li->fn_init = dlsym(li->lib_handle, "light_dll_init");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	li->fn_solve_lam_slice = dlsym(li->lib_handle, "light_dll_solve_lam_slice");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	li->light_ver = dlsym(li->lib_handle, "light_dll_ver");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}



(*li->light_ver)(sim);
(*li->fn_init)(sim);
}


void light_free_dlls(struct simulation *sim,struct light *li)
{
dlclose(li->lib_handle);
}

void light_free(struct simulation *sim,struct light *li)
{
	light_free_memory(sim,li);
	light_free_dlls(sim,li);
}
