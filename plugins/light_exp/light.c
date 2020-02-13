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
	@brief Exponential light solver.
*/


#include <util.h>
#include <dump_ctrl.h>
#include <gpvdm_const.h>
#include <light.h>
#include <device.h>
#include <light_interface.h>

#include <functions.h>
#include <log.h>


EXPORT void light_dll_ver(struct simulation *sim)
{
        printf_log(sim,"Exponential light model\n");
}

EXPORT int light_dll_solve_lam_slice(struct simulation *sim,struct device *dev,struct light *li, int z, int x, int l)
{
	int y;
	char temp[100];

	struct dim_light *dim=&li->dim;

	if (li->sun_E[l]==0.0)
	{
		return 0;
	}

	if (get_dump_status(sim,dump_optics)==TRUE)
	{
		sprintf(temp,"Solve light optical slice at %Lf nm\n",dim->l[l]*1e9);

		waveprint(sim,temp,dim->l[l]*1e9);
	}


	long double complex n0=0.0+0.0*I;

	//complex long double r=0.0+0.0*I;
	complex long double t=0.0+0.0*I;
	long double complex beta0=0.0+0.0*I;
	complex long double Ep=li->sun_E[l]+0.0*I;
	complex long double En=0.0+0.0*I;
	long double dy=dim->y[1]-dim->y[0];

	for (y=0;y<dim->ylen;y++)
	{

		n0=li->nbar[z][x][y][l];
		beta0=(2*PI*n0/dim->l[l]);
		Ep=Ep*cexp(-beta0*dy*I);

		t=li->t[z][x][y][l];



		li->Ep[z][x][y][l]=creal(Ep);
		li->Epz[z][x][y][l]=cimag(Ep);
		li->En[z][x][y][l]=creal(En);
		li->Enz[z][x][y][l]=cimag(En);

		if (y!=(dim->ylen-1))
		{
			if ((li->n[z][x][y][l]!=li->n[z][x][y+1][l])||(li->alpha[z][x][y][l]!=li->alpha[z][x][y+1][l]))
			{
				Ep=Ep*t;
			}
		}
	}

return 0;
}


