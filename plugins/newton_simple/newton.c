//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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


