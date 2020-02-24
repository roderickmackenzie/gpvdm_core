// 
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
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
// 

/** @file mesh.h
@brief meshing functions
*/

#ifndef mesh_h
#define mesh_h

void mesh_check_y(struct simulation *sim,struct mesh *in,struct device *dev);
void mesh_remesh(struct simulation *sim,struct mesh *in,struct device *dev);
void mesh_save(struct simulation *sim,char *file_name,struct mesh *in);
void mesh_free(struct mesh *in);
void mesh_build(struct simulation *sim,struct device *in);
void mesh_cal_layer_widths(struct device *in);
void mesh_init(struct mesh *in);
void mesh_load_file(struct simulation * sim, struct mesh *in,char *file);
void mesh_cpy(struct simulation *sim,struct mesh *out,struct mesh *in);
void mesh_malloc_sub_mesh(struct simulation * sim, struct mesh *in);
void mesh_gen_simple(struct simulation * sim, struct mesh *in,long double len,int points);
long double mesh_to_dim(struct simulation *sim,struct dimensions *dim, struct mesh *in,char xyz);
long double mesh_to_dim_heat(struct simulation *sim,struct dim_heat *dim, struct mesh *in,char xyz);
long double mesh_to_lin_array(struct simulation *sim,long double *mesh, long double *dmesh, struct mesh *in);

//mesh obj
void mesh_obj_cpy(struct simulation *sim,struct mesh_obj *out,struct mesh_obj *in);
void mesh_obj_load(struct simulation *sim,struct mesh_obj *mesh);
void mesh_obj_free(struct simulation *sim,struct mesh_obj *in);
void mesh_obj_init(struct mesh_obj *in);
void mesh_numerate_points(struct simulation *sim,struct device *in);
void mesh_dump(struct simulation *sim,struct dimensions *dim);
void mesh_dump_y(struct simulation *sim,struct dimensions *dim);
#endif
