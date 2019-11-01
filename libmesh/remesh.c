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

/** @file remesh.c
	@brief Remesh the solver.
*/


#include "sim.h"
#include "dump.h"
#include "remesh.h"
#include "newton_tricks.h"


struct remesh my_mesh;

void remesh_shrink(struct simulation *sim,struct device *in)
{

	if (in->dynamic_mesh==TRUE)
	{
		mesh_obj_cpy(sim,&(in->mesh_data_save),&(in->mesh_data));
		mesh_free(&(in->mesh_data.meshdata_x));
		mesh_gen_simple(sim, &(in->mesh_data.meshdata_x),in->xlen,4);
		mesh_to_dim(sim,&(in->dim), &(in->mesh_data.meshdata_x),'x');

	}

}

void remesh_expand_array_band(gdouble **y,int band,struct device *in)
{
}

void remesh_expand_array(gdouble *y,struct device *in)
{
}

void remesh_reset(struct simulation *sim,struct device *in,gdouble voltage)
{
}
