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
	@brief Flat optical profile solver.
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
        printf_log(sim,"Flat light model\n");
}

EXPORT int light_dll_solve_lam_slice(struct simulation *sim,struct device *dev,struct light *li, int z, int x, int l)
{
	int y;
	char one[100];
	struct dim_light *dim=&li->dim;
	struct epitaxy *epi=&(dev->my_epitaxy);

	if (li->sun_E[l]==0.0)
	{
		return 0;
	}

	if (get_dump_status(sim,dump_optics)==TRUE)
	{
		sprintf(one,"Solve light optical slice at %Lf nm\n",dim->l[l]*1e9);
		waveprint(sim,one,dim->l[l]*1e9);
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

		if ((dim->y[y]>epi->device_start) &&(dim->y[y]<epi->device_stop))
		{
			beta0=creal(beta0)+I*0;
		}

		Ep=Ep*cexp(-beta0*dy*I);

		//r=in->r[lam][i];
		t=li->t[z][x][y][l];

		//if ((in->n[lam][i]!=in->n[lam][i+1])||(in->alpha[lam][i]!=in->alpha[lam][i+1]))
		//{
		//	En=Ep*r;
		//}

		li->Ep[z][x][y][l]=creal(Ep);
		li->Epz[z][x][y][l]=cimag(Ep);
		li->En[z][x][y][l]=creal(En);
		li->Enz[z][x][y][l]=cimag(En);

		if ((li->n[z][x][y][l]!=li->n[z][x][y+1][l])||(li->alpha[z][x][y][l]!=li->alpha[z][x][y+1][l]))
		{
			Ep=Ep*t;
		}
	}

	return 0;
}


