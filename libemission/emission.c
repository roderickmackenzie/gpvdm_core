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

/** @file pl.c
	@brief Peform PL spectra.
*/


#include <stdio.h>
#include <dump.h>
#include <string.h>
#include <exp.h>
#include <dos.h>
#include "sim.h"
#include "i.h"
#include "dat_file.h"
#include "pl.h"
#include <cal_path.h>
#include <lang.h>
#include <color.h>



long double calculate_photon_power_m2(struct simulation *sim,struct device *in)
{

int y=0;
int layer=0;
long double lmax=0.0;
long double tot=0.0;
long double E=0.0;
long double R=0.0;
long double eff=0.0;
struct dimensions *dim=&(in->ns.dim);

	for (y=0;y<dim->ylen;y++)
	{
			layer=in->imat_epitaxy[0][0][y];
			if (in->my_epitaxy.layer[layer].pl_enabled==TRUE)
			{
				E=hp*cl/in->my_epitaxy.layer[layer].peak_wavelength;
				R=dim->dy[y]*(in->Rn[0][0][y]+in->Rp[0][0][y])/2.0;
				eff=in->my_epitaxy.layer[layer].avg_photon_extract_eff;
				tot+=R*E*eff;

				//printf("%Le %Le %Le %Le %Le\n",R,E,eff,in->Rn[0][0][y]);
				//getchar();

			}
	}

	//getchar();
return tot;
}



