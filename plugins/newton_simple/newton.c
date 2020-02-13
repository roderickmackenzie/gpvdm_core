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
#include <gpvdm_const.h>
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
#include <circuit.h>




void fill_matrix(struct simulation *sim,struct device *in,int z,int x_in)
{
int i=0;
int x=0;
int x_max=0;
int band=0;
struct dimensions *dim=&in->ns.dim;
struct circuit *cir=&(in->cir);


		/*
		long double n0=1.0;
		long double J0=1.0;

		long double J=0.0;
		long double V=0.0;
		long double J_light=0.0;

		J_light=0.0;

		for (i=0;i<dim->ylen;i++)
		{
			J_light+=(in->Gn[z][x][i]+in->Gp[z][x][i])*dim->dy[i]/2.0;
		}

		V=in->Vapplied_y0[z][x]-in->Vapplied_y1[z][x];
		J=(J0/Q)*(expl(((Q*V)/(n0*kb*in->Tl[z][x][0])))-1.0)-J_light;
		*/

	circuit_solve(sim,cir,in);
	circuit_transfer_to_electrical_mesh(sim,cir,in);

}

void dllinternal_solver_realloc(struct simulation *sim,struct device *in,int dim)
{
}

void dllinternal_solver_free_memory(struct simulation *sim,struct device *in)
{
}

int dllinternal_solve_cur(struct simulation *sim,struct device *in, int z, int x)
{

	fill_matrix(sim,in,z,x);

	struct newton_state *ns=&(in->ns);

	ns->last_error=0.0;//error;
	ns->last_ittr=0.0;

	in->newton_last_ittr=0;

	in->odes+=1;


return 0;
}


