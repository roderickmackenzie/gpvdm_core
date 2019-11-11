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

/** @file light.c
	@brief Read an optical generation profile from a file.
*/


#include <util.h>
#include <dump_ctrl.h>
#include <complex_solver.h>
#include <const.h>
#include <light.h>
#include <device.h>
#include <light_interface.h>
#include <string.h>
#include <functions.h>
#include <log.h>
#include <inp.h>
#include <cal_path.h>

EXPORT void light_dll_ver(struct simulation *sim)
{
        printf_log(sim,"External field light model\n");
}

EXPORT int light_dll_solve_lam_slice(struct simulation *sim,struct device *cell,struct light *in,int lam)
{
	int i=0;
	long double G=0.0;
	struct inp_file inp;
	char temp[PATH_MAX];

	inp_init(sim,&inp);

	inp_load_from_path(sim,&inp,get_input_path(sim),"light.inp");

	inp_search_gdouble(sim,&inp,&(G),"#light_flat_generation_rate");

	inp_free(sim,&inp);

	sprintf(temp,"Using flat generation rate of %Le\n",G);



	G=G*in->Psun;

	in->disable_cal_photon_density=TRUE;

	for (i=0;i<in->points;i++)
	{
		in->Gn[i]=G;
		in->Gp[i]=G;
	}

	in->finished_solveing=TRUE;

return 0;
}


