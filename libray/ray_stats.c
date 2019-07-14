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

#include <stdio.h>
#include <ray.h>
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <ray_fun.h>

/** @file ray_stats.c
	@brief Perfrom stats on the ray tracing image
*/
void ray_escape_angle_reset(struct image *in,int l)
{
	int i;

	for (i=0;i<in->escape_bins;i++)
	{
		in->ang_escape[l][i]=0.0;
	}
}

void ray_cal_escape_angle(struct image *in,int l)
{
	int i;
	int ii;
	double mag=0.0;

	for (i=0;i<in->nrays;i++)
	{
		if (in->rays[i].state==DONE)
		{
			if (in->rays[i].xy_end.y<in->y_escape_level)
			{
				mag=in->rays[i].mag;
				double raw_ang=(360.0/(2*3.14159))*atan(fabs(in->rays[i].dir.y)/in->rays[i].dir.x);
				if (raw_ang<0.0)
				{
					raw_ang=raw_ang+180.0;
				}
				double ang=raw_ang;
				int bin=(int)((ang/180.0)*(double)in->escape_bins);
				in->ang_escape[l][bin]+=in->rays[i].mag;
				//printf("%lf %lf %lf %lf\n",in->rays[i].dir.x,fabs(in->rays[i].dir.y),ang,raw_ang);
				//getchar();
			}
		}
		
	}

}

double get_eff(struct image *in)
{
int i;
double tot=0.0;
	for (i=0;i<in->nrays;i++)
	{
		if (in->rays[i].state==DONE)
		{
			if (in->rays[i].xy_end.y<in->y_escape_level)
			{
				tot+=in->rays[i].mag;
			}
		}
		
	}

return tot;
}
