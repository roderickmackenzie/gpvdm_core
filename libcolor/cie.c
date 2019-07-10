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
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <color.h>

/** @file wavelength_to_rgb.c
	@brief Turn a wavelegngth to an RGB color.
*/

void color_cie_load(struct simulation *sim,struct device *in)
{
	char path[PATH_MAX];
	join_path(2, path,get_cie_color_path(sim),"x.inp");

	inter_load(sim,&(in->cie_x),path);
	inter_sort(&(in->cie_x));

	join_path(2, path,get_cie_color_path(sim),"y.inp");
	inter_load(sim,&(in->cie_y),path);
	inter_sort(&(in->cie_y));

	join_path(2, path,get_cie_color_path(sim),"z.inp");
	inter_load(sim,&(in->cie_z),path);
	inter_sort(&(in->cie_z));
}

void color_cie_cal_XYZ_from_eV(struct simulation *sim,struct device *in,long double *X,long double *Y,long double *Z,struct istruct *L_eV)
{
int i;
struct istruct L;
inter_copy(&L,L_eV,TRUE);
*X=0.0;
*Y=0.0;
*Z=0.0;

for (i=0;i<L.len;i++)
{
	L.x[i]=	(hp*cl)/(L.x[i]*Q);
}

inter_sort(&L);

inter_save(&L,"e.dat");

long double dl=0.0;
//inter_dump(sim,&(in->cie_x));
for (i=0;i<in->cie_x.len-1;i++)
{
	dl=in->cie_x.x[i+1]-in->cie_x.x[i];
	(*X)+=dl*inter_get_hard(&L,in->cie_x.x[i])*(in->cie_x.data[i]);

	dl=in->cie_y.x[i+1]-in->cie_y.x[i];
	(*Y)+=dl*inter_get_hard(&L,in->cie_y.x[i])*(in->cie_y.data[i]);

	dl=in->cie_z.x[i+1]-in->cie_z.x[i];
	(*Z)+=dl*inter_get_hard(&L,in->cie_z.x[i])*(in->cie_z.data[i]);
	//printf("X=%Le %Le %Le\n",inter_get_hard(&L,in->cie_y.x[i]),in->cie_x.data[i],dl);
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
	r=255.0;
	g=0.0;
	b=0.0;

}
//printf("%Le %Le %Le\n",r,g,b);
*R=(int)(r);
*G=(int)(g);
*B=(int)(b);
}

void color_cie_free(struct simulation *sim,struct device *in)
{
	inter_free(&(in->cie_x));
	inter_free(&(in->cie_y));
	inter_free(&(in->cie_z));

}
