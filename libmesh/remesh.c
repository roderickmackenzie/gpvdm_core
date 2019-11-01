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
#include <mesh.h>
#include <memory.h>

struct remesh my_mesh;

void restore_state(struct simulation *sim,struct device *in)
{
	int band=0;
	int x=0;
	int z=0;

	three_d_interpolate_gdouble(in->ns.phi, in->ns_save.phi, &(in->dim),&(in->dim_save));

	three_d_interpolate_gdouble(in->ns.x, in->ns_save.x, &(in->dim),&(in->dim_save));
	three_d_interpolate_gdouble(in->ns.xp, in->ns_save.xp, &(in->dim),&(in->dim_save));

	for (band=0;band<in->dim.srh_bands;band++)
	{
		three_d_interpolate_srh(in->ns.xt, in->ns_save.xt, &(in->dim),&(in->dim_save),band);
		three_d_interpolate_srh(in->ns.xpt, in->ns_save.xpt, &(in->dim),&(in->dim_save),band);
	}


	for (z=0;z<in->dim.zmeshpoints;z++)
	{
		for (x=0;x<in->dim.xmeshpoints;x++)
		{
			update_y_array(sim,in,z,x);
		}
	}
}

void remesh_shrink(struct simulation *sim,struct device *in)
{

	if (in->dynamic_mesh==TRUE)
	{
		dim_cpy(&(in->dim_save),&(in->dim));
		newton_save_state_cpy(&(in->ns_save),&(in->ns),&(in->dim));
		mesh_obj_cpy(sim,&(in->mesh_data_save),&(in->mesh_data));
		mesh_free(&(in->mesh_data.meshdata_x));
		mesh_gen_simple(sim, &(in->mesh_data.meshdata_x),in->xlen,10);
		mesh_to_dim(sim,&(in->dim), &(in->mesh_data.meshdata_x),'x');

		restore_state(sim,in);

	}

}

void remesh_expand_array_band(gdouble **y,int band,struct device *in)
{
}


void remesh_reset(struct simulation *sim,struct device *in,gdouble voltage)
{

	if (in->dynamic_mesh==TRUE)
	{
		newton_save_state_cpy(&(in->ns_save),&(in->ns),&(in->dim));
		dim_swap(&(in->dim),&(in->dim_save));
		mesh_obj_cpy(sim,&(in->mesh_data_save),&(in->mesh_data));

		restore_state(sim,in);
		sim_externalv(sim,in,voltage);
		//three_d_quick_dump("0.dat", in->ns.phi, &(in->dim));
		//three_d_quick_dump("1.dat", in->ns_save.phi, &(in->dim_save));
		//getchar();
		/*
		newton_push_state(in);
		in->max_electrical_itt=1000;
		in->newton_min_itt=10;
		in->electrical_clamp=4.0;

		//save the old newton state (done)
		//save the old dim (done)
		//restore the mesh (done)
		//interpolate result back to



		for (band=0;band<in->srh_bands;band++)
		{
			remesh_expand_array_band(in->xt,band,in);
			remesh_expand_array_band(in->xpt,band,in);
		}


		for (i=0;i<my_mesh.len;i++)
		{
			in->ymesh[i]=my_mesh.x[i];
		}


		in->ymeshpoints=my_mesh.len;
		my_mesh.len=0;
		update_arrays(in);

		sim_externalv(in,voltage);
		newton_pop_state(in);

		free(my_mesh.x);
		*/
	}
}
