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
#include <gpvdm_const.h>
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

EXPORT int light_dll_solve_lam_slice(struct simulation *sim,struct device *cell,struct light *li,int z, int x, int l)
{
	int y=0;
	struct dim_light *dim=&li->dim;
	long double G=0.0;

	G=G*li->Psun;

	li->disable_cal_photon_density=TRUE;

	for (y=0;y<dim->ylen;y++)
	{
		li->Gn[z][x][y]=li->light_flat_generation_rate;
		li->Gp[z][x][y]=li->light_flat_generation_rate;
	}

	li->finished_solveing=TRUE;

return 0;
}


