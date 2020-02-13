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
#include <i.h>
#include <gpvdm_const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <color.h>

/** @file wavelength_to_rgb.c
	@brief Turn a wavelegngth to an RGB color.
*/

void color_cie_init(struct simulation *sim)
{
	inter_init(sim,&(sim->cie_x));
	inter_init(sim,&(sim->cie_y));
	inter_init(sim,&(sim->cie_z));

}

void color_cie_load(struct simulation *sim)
{
	char path[PATH_MAX];
	join_path(2, path,get_cie_color_path(sim),"x.inp");

	inter_load(sim,&(sim->cie_x),path);
	inter_sort(&(sim->cie_x));

	join_path(2, path,get_cie_color_path(sim),"y.inp");
	inter_load(sim,&(sim->cie_y),path);
	inter_sort(&(sim->cie_y));

	join_path(2, path,get_cie_color_path(sim),"z.inp");
	inter_load(sim,&(sim->cie_z),path);
	inter_sort(&(sim->cie_z));

}

void color_cie_cal_XYZ(struct simulation *sim,long double *X,long double *Y,long double *Z,struct istruct *L_input, int input_in_ev)
{
int i;
struct istruct L;
inter_init(sim,&L);
inter_copy(&L,L_input,TRUE);
*X=0.0;
*Y=0.0;
*Z=0.0;

if (input_in_ev==TRUE)
{
	for (i=0;i<L.len;i++)
	{
		L.x[i]=	(hp*cl)/(L.x[i]*Q);
	}
}
inter_sort(&L);

//inter_save(&L,"e.dat");

long double dl=0.0;
//printf("urhu\n");
//inter_dump(sim,&(sim->cie_x));

if (sim->cie_x.len==0)
{
	printf("sim->cie_x.len==0\n");
	exit(0);
}


for (i=0;i<sim->cie_x.len-1;i++)
{

	dl=sim->cie_x.x[i+1]-sim->cie_x.x[i];
	(*X)+=dl*inter_get_hard(&L,sim->cie_x.x[i])*(sim->cie_x.data[i]);

	dl=sim->cie_y.x[i+1]-sim->cie_y.x[i];
	(*Y)+=dl*inter_get_hard(&L,sim->cie_y.x[i])*(sim->cie_y.data[i]);

	dl=sim->cie_z.x[i+1]-sim->cie_z.x[i];
	(*Z)+=dl*inter_get_hard(&L,sim->cie_z.x[i])*(sim->cie_z.data[i]);
	//printf("X=%Le %Le %Le\n",inter_get_hard(&L,sim->cie_y.x[i]),sim->cie_x.data[i],dl);
	//getchar();
}

//printf("SUM>>>%Le\n",*X);
//printf("SUM>>>%Le\n",*Y);
//printf("SUM>>>%Le\n",*Z);
inter_free(&L);
}

void color_XYZ_to_rgb(int *R,int *G, int *B,long double X,long double Y, long double Z)
{
long double r=0.0;
long double g=0.0;
long double b=0.0;
long double max=0.0;

r=0.41847*X-0.15866*Y-0.082835*Z;
g=-0.091169*X+0.25243*Y+0.015708*Z;
b=0.00092090*X-0.0025498*Y+0.17860*Z;

if (r<0.0) r=0.0;
if (g<0.0) g=0.0;
if (b<0.0)  b=0.0;

max=0.0;
if (r>max) max=r;
if (g>max) max=g;
if (b>max) max=b;

if (max!=0)
{
	r=255.0*r/max;
	g=255.0*g/max;
	b=255.0*b/max;
}else
{
	r=0.0;
	g=0.0;
	b=0.0;

}
//printf("rgb  %Le %Le %Le\n",r,g,b);
*R=(int)(r);
*G=(int)(g);
*B=(int)(b);
}

void color_XYZ_to_rgb_no_norm(long double *R,long double *G, long double *B,long double X,long double Y, long double Z)
{
long double r=0.0;
long double g=0.0;
long double b=0.0;
long double max=0.0;

r=0.41847*X-0.15866*Y-0.082835*Z;
g=-0.091169*X+0.25243*Y+0.015708*Z;
b=0.00092090*X-0.0025498*Y+0.17860*Z;

if (r<0.0) r=0.0;
if (g<0.0) g=0.0;
if (b<0.0)  b=0.0;

max=0.0;
if (r>max) max=r;
if (g>max) max=g;
if (b>max) max=b;

if (max!=0)
{
	r=255.0*r/max;
	g=255.0*g/max;
	b=255.0*b/max;
}else
{
	r=0.0;
	g=0.0;
	b=0.0;

}
//printf("rgb  %Le %Le %Le\n",r,g,b);
*R=(int)(r);
*G=(int)(g);
*B=(int)(b);
}
void color_cie_free(struct simulation *sim)
{
	inter_free(&(sim->cie_x));
	inter_free(&(sim->cie_y));
	inter_free(&(sim->cie_z));

}

