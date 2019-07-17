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

/** @file ray_shapes.c
	@brief Basic shapes for ray tracing
*/

void add_triangle(struct image *in, double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,int id,int edge)
{
	in->p[in->lines].xy0.x=x0;
	in->p[in->lines].xy0.y=y0;
	in->p[in->lines].xy0.z=z0;

	in->p[in->lines].xy1.x=x1;
	in->p[in->lines].xy1.y=y1;
	in->p[in->lines].xy1.z=z1;

	in->p[in->lines].xy2.x=x2;
	in->p[in->lines].xy2.y=y2;
	in->p[in->lines].xy2.z=z2;

	in->p[in->lines].edge=edge;
	in->p[in->lines].id=id;
	in->lines++;
}


void add_box(struct image *in,double x0,double y0,double z0,double dx,double dy,double dz,int material,int sim_edge)
{
	//btm
	add_triangle(in,
					x0,y0,z0,
					x0+dx,y0,
					z0,x0,y0,z0+dz,
													in->objects,sim_edge);
	add_triangle(in,
					x0+dx	,	y0,z0,
					x0+dx	,	y0,z0+dz,
					x0   	,	y0,z0+dz,
													in->objects,sim_edge);

	//top
	add_triangle(in,
					x0		,y0+dy	,	z0		,
					x0+dx	,y0+dy	,	z0		,
					x0		,y0+dy	,	z0+dz	,
													in->objects,sim_edge);
	add_triangle(in,
					x0+dx	,y0+dy	,z0			,
					x0+dx	,y0+dy	,z0+dz		,
					x0   	,y0+dy	,z0+dz		,
													in->objects,sim_edge);
	
	sim_edge=TRUE;
	//left
	add_triangle(in,
					x0		,y0		,z0			,
					x0		,y0+dy	,z0			,
					x0   	,y0		,z0+dz		,
													in->objects,sim_edge);

	add_triangle(in,
					x0		,y0+dy	,z0			,
					x0		,y0+dy	,z0+dz		,
					x0   	,y0		,z0+dz		,
													in->objects,sim_edge);

	//right
	add_triangle(in,
					x0+dx		,y0		,z0			,
					x0+dx		,y0+dy	,z0			,
					x0+dx   	,y0		,z0+dz		,
													in->objects,sim_edge);

	add_triangle(in,
					x0+dx	,y0+dy	,z0			,
					x0+dx	,y0+dy	,z0+dz		,
					x0+dx  	,y0		,z0+dz		,
													in->objects,sim_edge);

	//front
	add_triangle(in,
					x0			,y0		,z0		,
					x0+dx		,y0		,z0		,
					x0   		,y0+dy	,z0		,
													in->objects,sim_edge);

	add_triangle(in,
					x0			,y0+dy	,z0		,
					x0+dx		,y0+dy	,z0		,
					x0+dx   	,y0		,z0		,
													in->objects,sim_edge);
					
	//back
	add_triangle(in,
					x0			,y0		,z0+dz		,
					x0+dx		,y0		,z0+dz		,
					x0   		,y0+dy	,z0+dz		,
													in->objects,sim_edge);

	add_triangle(in,
					x0			,y0+dy	,z0+dz		,
					x0+dx		,y0+dy	,z0+dz		,
					x0+dx   	,y0		,z0+dz		,
													in->objects,sim_edge);
	
	//top
	in->obj_mat_number[in->objects]=material;
	in->objects++;
}

