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
#include "const.h"
#include "device.h"
#include "dump_ctrl.h"
#include "config.h"
#include "cal_path.h"
#include "lang.h"
#include "log.h"
#include "sim.h"
#include "memory.h"

static int unused __attribute__((unused));

static gdouble last_Psun= -1000.0;
static gdouble last_laser_eff= -1000.0;
static gdouble last_wavelength_laser= -1000.0;

void light_load_dlls(struct simulation *sim,struct light *in)
{
	char lib_path[200];
	char lib_name[100];

	printf_log(sim,"%s\n",_("Initializing optical model"));

	sprintf(lib_name,"light_%s",in->mode);
	find_dll(sim, lib_path,lib_name);


	char *error;

	in->lib_handle = dlopen(lib_path, RTLD_LAZY);

	if (!in->lib_handle)
	{
		ewe(sim, "%s\n", dlerror());
	}

	in->fn_init = dlsym(in->lib_handle, "light_dll_init");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	in->fn_solve_lam_slice = dlsym(in->lib_handle, "light_dll_solve_lam_slice");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}

	in->light_ver = dlsym(in->lib_handle, "light_dll_ver");
	if ((error = dlerror()) != NULL)
	{
		ewe(sim, "%s\n", error);
	}



(*in->light_ver)(sim);
(*in->fn_init)(sim);
}

void light_solve_and_update(struct simulation *sim,struct device *cell,struct light *in,gdouble laser_eff_in)
{
	int x=0;
	int y=0;
	int z=0;

	if (in->disable_transfer_to_electrical_mesh==FALSE)
	{
		if (fabs(in->device_ylen-cell->ylen)>0.01e-9)
		{
		ewe(sim,"The electrical mesh (%.9le) and the optical mesh (%.9le) don't match. %le",cell->ylen,in->device_ylen);
		}
	}


	in->laser_eff=laser_eff_in;
	if ((last_laser_eff!=in->laser_eff)||(last_Psun!=in->Psun)||(last_wavelength_laser!=in->laser_wavelength)||(in->force_update==TRUE))
	{
		light_solve_optical_problem(sim,in);
		last_laser_eff=in->laser_eff;
		last_Psun=in->Psun;
		last_wavelength_laser=in->laser_wavelength;
	}

	if (in->laser_pos!=-1)
	{
		light_dump_1d(sim,in, in->laser_pos,"");
	}

	light_transfer_gen_rate_to_device(cell,in);

	if (in->flip_field==TRUE)
	{
		gdouble ***Gn;
		gdouble ***Gp;

		malloc_3d_gdouble(cell, &Gn);
		malloc_3d_gdouble(cell, &Gp);

		for (z=0;z<cell->zmeshpoints;z++)
		{

			for (x=0;x<cell->xmeshpoints;x++)
			{

				for (y=0;y<cell->ymeshpoints;y++)
				{
					Gn[z][x][y]=cell->Gn[z][x][y];
					Gp[z][x][y]=cell->Gp[z][x][y];
				}

			}
		}

		for (z=0;z<cell->zmeshpoints;z++)
		{
			for (x=0;x<cell->xmeshpoints;x++)
			{
				for (y=0;y<cell->ymeshpoints;y++)
				{
					cell->Gn[z][x][cell->ymeshpoints-y-1]=Gn[z][x][y];
					cell->Gp[z][x][cell->ymeshpoints-y-1]=Gp[z][x][y];
				}
			}
		}

		for (z=0;z<cell->zmeshpoints;z++)
		{
			for (x=0;x<cell->xmeshpoints;x++)
			{
				for (y=0;y<cell->ymeshpoints;y++)
				{
					cell->Gn[z][x][cell->ymeshpoints-y-1]=Gn[z][x][y];
					cell->Gp[z][x][cell->ymeshpoints-y-1]=Gp[z][x][y];
				}
			}
		}

		free_3d_gdouble(cell, Gn);
		free_3d_gdouble(cell, Gp);
	}

}

void light_init(struct light *in)
{
	last_Psun= -1.0;
	last_laser_eff= -1.0;
	last_wavelength_laser= -1.0;
	in->laser_wavelength= -1.0;
	in->laser_pos= -1;
	in->laser_wavelength= -1.0;
	in->lstart= -1.0;
	in->spotx= -1.0;
	in->spoty= -1.0;
	in->pulseJ= -1.0;
	in->pulse_width= -1.0;
	in->layers=-1;
	in->G_percent=NULL;
	in->disable_cal_photon_density=FALSE;
}


void light_load_config(struct simulation *sim,struct light *in,struct epitaxy *my_epitaxy)
{
	light_load_config_file(sim,in);
	light_import_epitaxy(sim,in,my_epitaxy);
	light_load_materials(sim,in);
	light_memory(sim,in);
	light_init_mesh(sim,in,my_epitaxy);
}

int light_solve_lam_slice(struct simulation *sim, struct light *in,int lam)
{
in->disable_cal_photon_density=FALSE;
return (*in->fn_solve_lam_slice)(sim,in,lam);
}


void light_free_dlls(struct simulation *sim,struct light *in)
{
dlclose(in->lib_handle);
}

void light_free(struct simulation *sim,struct light *in)
{
light_free_memory(sim,in);
light_free_dlls(sim,in);
}
