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

/** @file mesh_struct.h
@brief meshing structure
*/

#ifndef mesh_struct_h
#define mesh_struct_h

struct mesh_layer
{
	long double dx;
	long double len;
	long double mul;
	long double *dmesh;
	int n_points;
	int left_right;
};

struct mesh
{
	long double start;
	struct mesh_layer *layers;
	int nlayers;
	int remesh;
	int tot_points;
};

struct mesh_obj
{
	struct mesh meshdata_x;
	struct mesh meshdata_y;
	struct mesh meshdata_z;
};

#endif
