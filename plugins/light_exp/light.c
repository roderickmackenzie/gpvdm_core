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


