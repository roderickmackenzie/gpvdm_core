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

/** @file light_utils.c
	@brief Miscellaneous functions for the light model.
*/

#include "util.h"
#include "gpvdm_const.h"
#include "light.h"
#include "device.h"
#include "gpvdm_const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "lang.h"
#include "log.h"
#include "memory.h"
#include <light_fun.h>

static int unused __attribute__((unused));

void light_solve_and_update(struct simulation *sim,struct device *dev,struct light *li,long double laser_eff_in)
{

	struct dimensions *dim=&dev->ns.dim;

	li->laser_eff=laser_eff_in;

	if ((li->last_laser_eff!=li->laser_eff)||(li->last_Psun!=li->Psun)||(li->last_wavelength_laser!=li->laser_wavelength)||(li->force_update==TRUE))
	{
		light_solve_optical_problem(sim,dev,li);
		li->last_laser_eff=li->laser_eff;
		li->last_Psun=li->Psun;
		li->last_wavelength_laser=li->laser_wavelength;
	}

	if (li->laser_pos!=-1)
	{
		light_dump_1d(sim,li, li->laser_pos,"");
	}

	light_transfer_gen_rate_to_device(dev,li);

	if (li->flip_field==TRUE)
	{
		flip_zxy_long_double_y(sim, dim,dev->Gn);
		flip_zxy_long_double_y(sim, dim,dev->Gp);
	}

}

void light_solve_optical_problem(struct simulation *sim,struct device *dev,struct light *li)
{
	int l=0;
	long double Psun=0.0;
	struct dim_light *dim=&li->dim;

	Psun=li->Psun*gpow(10.0,-li->ND);

	light_set_sun_power(li,Psun,li->laser_eff);

	light_solve_all(sim,dev,li);

	light_cal_photon_density(sim,li,dev);

	for (l=0;l<dim->llen;l++)
	{
		light_dump_1d(sim,li, l,"");
	}

}

void light_solve_all(struct simulation *sim,struct device *cell,struct light *li)
{
	int z=0;
	int x=0;
	int l=0;

	int slices_solved=0;
	struct dim_light *dim=&li->dim;
	li->finished_solveing=FALSE;
	char temp[200];

	for (l=0;l<dim->llen;l++)
	{
		for (z=0;z<dim->zlen;z++)
		{
			for (x=0;x<dim->xlen;x++)
			{

				//printf("here %Le %d\n",li->sun_E[l],l);
				memset_light_zxyl_long_double_y(dim, li->En,z,x,l,0.0);
				memset_light_zxyl_long_double_y(dim, li->Ep,z,x,l,0.0);
				memset_light_zxyl_long_double_y(dim, li->Enz,z,x,l,0.0);
				memset_light_zxyl_long_double_y(dim, li->Epz,z,x,l,0.0);

				if (li->sun_E[l]!=0.0)
				{
					if ((x==0)&&(z==0))
					{
						if (li->print_wavlengths==TRUE)
						{
							if (get_dump_status(sim,dump_optics)==TRUE)
							{
								sprintf(temp,"%s %Lf nm\n",_("Solve optical slice at"),dim->l[l]*1e9);
								waveprint(sim,temp,dim->l[l]*1e9);
							}
						}
					}

					light_solve_lam_slice(sim,cell,li,z,x,l);

					//printf("here %d\n",li->finished_solveing);
					//getchar();

					if (li->finished_solveing==TRUE)
					{
						break;
					}

					slices_solved++;
				}//else
				//{

				//}
			}
		}
	}

	if (slices_solved==0)
	{
		printf_log(sim,_("It was dark I did not solve any slices\n"));
	}

}


int light_solve_lam_slice(struct simulation *sim, struct device *dev, struct light *li,int z, int x,int l)
{
li->disable_cal_photon_density=FALSE;
return (*li->fn_solve_lam_slice)(sim,dev,li,z,x,l);
}
