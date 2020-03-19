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


/** @file device_heat.c
	@brief Malloc for the device structure.
*/

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <device.h>
#include <string.h>
#include <dump.h>
#include <mesh.h>
#include <ray_fun.h>
#include <newton_tricks.h>
#include <memory.h>
#include <circuit.h>
#include <shape.h>
#include <lang.h>
#include <util.h>
#include <heat_fun.h>
#include <gpvdm_const.h>

static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));

void device_calculate_joule_heat(struct simulation *sim,struct device *dev)
{
	int z=0;
	int x=0;
	int y=0;
	long double Ecl=0.0;
	long double Ecr=0.0;
	long double Evl=0.0;
	long double Evr=0.0;
	long double yl=0.0;
	long double yr=0.0;
	long double dh=0.0;
	long double Jn=0.0;
	long double Jp=0.0;
	long double H_joule_e=0.0;
	long double H_joule_h=0.0;

	struct heat *thermal=&(dev->thermal);
	struct dimensions *dim=&(dev->ns.dim);

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				if (thermal->joule_heating==TRUE)
				{
					if (y==0)
					{
						Ecl=-dev->Xi[z][x][0]-dev->Vapplied_y0[z][x];
						Evl=-dev->Xi[z][x][0]-dev->Vapplied_y0[z][x]-dev->Eg[z][x][0];
						yl=dim->ymesh[0]-(dim->ymesh[1]-dim->ymesh[0]);
					}else
					{
						Ecl=dev->Ec[z][x][y-1];
						Evl=dev->Ev[z][x][y-1];
						yl=dim->ymesh[y-1];
					}

					if (y==dim->ylen-1)
					{
						Ecr=-dev->Xi[z][x][y]-dev->V_y1[z][x];
						Evr=-dev->Xi[z][x][y]-dev->V_y1[z][x]-dev->Eg[z][x][y];
						yr=dim->ymesh[y]+(dim->ymesh[y]-dim->ymesh[y-1]);

					}else
					{
						Ecr=dev->Ec[z][x][y];
						Evr=dev->Ev[z][x][y];
						yr=dim->ymesh[y];
					}

					dh=yr-yl;
					Jn=dev->Jn[z][x][y];
					Jp=dev->Jp[z][x][y];

					H_joule_e=(Ecr-Ecl)*Jn/dh;
					H_joule_h=(Evr-Evl)*Jp/dh;

					dev->He[z][x][y]=H_joule_e;//+dev->Hasorb[i]/2.0;
					dev->Hh[z][x][y]=H_joule_h;//+dev->Hasorb[i]/2.0;


					dev->H_joule[z][x][y]=H_joule_e+H_joule_h;
					//printf("%ld %ld %ld %Le %Le %Le\n",z,x,y,H_joule_e,H_joule_h,Jn,Jp);
				}else
				{
					dev->H_joule[z][x][y]=0.0;
				}
				//printf(">%Le %Le\n",dim->ymesh[y],dev->H_joule[z][x][y]);
			}
		}
	}
	//getchar();
}


void device_calculate_recombination_heat(struct simulation *sim,struct device *dev)
{
	int z=0;
	int x=0;
	int y=0;
	long double Ecl=0.0;
	long double Ecr=0.0;
	long double Evl=0.0;
	long double Evr=0.0;
	long double yl=0.0;
	long double yr=0.0;
	long double dh=0.0;
	long double Jn=0.0;
	long double Jp=0.0;
	long double H_joule_e=0.0;
	long double H_joule_h=0.0;

	struct heat *thermal=&(dev->thermal);
	struct dimensions *dim=&(dev->ns.dim);

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				if (thermal->recombination_heating==TRUE)
				{
					dev->H_recombination[z][x][y]=dev->Eg[z][x][y]*(dev->Rn_srh[z][x][y]+dev->Rp_srh[z][x][y]+dev->Rfree[z][x][y])*Q;
				}else
				{
					dev->H_recombination[z][x][y]=0.0;
				}
			}
		}
	}
}
