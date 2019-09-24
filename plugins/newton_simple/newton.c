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

/** @file newton.c
	@brief Full 2D newton solver.
*/



#include <string.h>
#include <log.h>
#include <const.h>
#include "newton.h"
#include <dll_export.h>
#include <util.h>
#include <exp.h>
#include <advmath.h>
#include <dump.h>
#include <cal_path.h>
#include <dos.h>
#include <sim.h>
#include <solver_interface.h>
#include <contacts.h>





void fill_matrix(struct simulation *sim,struct device *in,int z,int x_in)
{
int i=0;
int x=0;
int x_max=0;
int band=0;

	if (x_in==-1)
	{
		x=0;
		x_max=in->xmeshpoints;
	}else
	{
		x=x_in;
		x_max=x_in;
	}

//contacts_update(sim,in);
//contacts_dump(sim,in);

long double n0=1.0;
long double J0=1.0;

long double J=0.0;
long double V=0.0;
long double J_light=0.0;


	for (i=0;i<in->my_epitaxy.layers;i++)
	{
		if (in->my_epitaxy.layer[i].electrical_layer==TRUE)
		{
			n0=in->my_epitaxy.layer[i].n;
			J0=in->my_epitaxy.layer[i].J0;
			break;
		}
	}



do
{

	J_light=0.0;

	for (i=0;i<in->ymeshpoints;i++)
	{
		J_light+=(in->Gn[z][x][i]+in->Gp[z][x][i])*in->dymesh[i]/2.0;
	}

	V=in->Vapplied_l[z][x]-in->Vapplied_r[z][x];
	J=(J0/Q)*(expl(((Q*V)/(n0*kb*in->Tl[z][x][0])))-1.0)-J_light;

	//printf("%Le %Le %Le %Le\n",n0,J0,J,V);

	
	for (i=0;i<in->ymeshpoints;i++)
	{


		if (i==0)
		{
			in->Jnleft[z][x]=J;
			in->Jpleft[z][x]=J;
		}

		if (i==in->ymeshpoints-1)
		{
			in->Jnright[z][x]=J;
			in->Jpright[z][x]=J;
		}

		in->Jn[z][x][i]=J;
		in->Jp[z][x][i]=J;

		in->Jn_x[z][x][i]=J;
		in->Jp_x[z][x][i]=J;


		in->nrelax[z][x][i]=0.0;
		in->ntrap_to_p[z][x][i]=0.0;
		in->prelax[z][x][i]=0.0;
		in->ptrap_to_n[z][x][i]=0.0;


		for (band=0;band<in->srh_bands;band++)
		{

			in->nt_r1[z][x][i][band]=0.0;
			in->nt_r2[z][x][i][band]=0.0;
			in->nt_r3[z][x][i][band]=0.0;
			in->nt_r4[z][x][i][band]=0.0;

			in->pt_r1[z][x][i][band]=0.0;
			in->pt_r2[z][x][i][band]=0.0;
			in->pt_r3[z][x][i][band]=0.0;
			in->pt_r4[z][x][i][band]=0.0;

		}

		in->Rn[z][x][i]=0.0;
		in->Rp[z][x][i]=0.0;

		in->Rn_srh[z][x][i]=0.0;
		in->Rp_srh[z][x][i]=0.0;

	}

	x++;

}while(x<x_max);



}

void dllinternal_solver_realloc(struct simulation *sim,struct device *in,int dim)
{
}

void dllinternal_solver_free_memory(struct device *in)
{
}

int dllinternal_solve_cur(struct simulation *sim,struct device *in, int z, int x)
{

	fill_matrix(sim,in,z,x);

	in->last_error=0.0;//error;
	in->last_ittr=0.0;

	in->newton_last_ittr=0;

	in->odes+=1;


return 0;
}


