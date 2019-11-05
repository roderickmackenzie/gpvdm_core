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

	for (y=0;y<dim->ymeshpoints;y++)
	{
			layer=in->imat_epitaxy[0][0][y];
			if (in->my_epitaxy.layer[layer].pl_enabled==TRUE)
			{
				E=hp*cl/in->my_epitaxy.layer[layer].peak_wavelength;
				R=dim->dymesh[y]*(in->Rn[0][0][y]+in->Rp[0][0][y])/2.0;
				eff=in->my_epitaxy.layer[layer].avg_photon_extract_eff;
				tot+=R*E*eff;
				//printf("%d %Le\n",layer,tot);
			}
	}

	//getchar();
return tot;
}



