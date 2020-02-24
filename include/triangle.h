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

/** @file triangle.h
@brief triangle structure
*/
#ifndef triangle_h
#define triangle_h

#include <vec.h>

struct triangle
{
	struct vec xy0;
	struct vec xy1;
	struct vec xy2;
	int object_uid;

	int object_type;

	//pre calculation
	struct vec edge1;
	struct vec edge2;

	int obj_left;
	int obj_right;

	//not used
	int tri_uid;		
};

struct triangles
{
	struct triangle *data;
	int edges_calculated;
	int max_len;
	int len;
};

void triangle_print(struct triangle *in);
void triangle_cog(struct vec *out,struct triangle *in);
void triangle_norm(struct vec *ret,struct triangle *my_obj);
void triangle_dump(char *file_name,struct triangle *tri);
int triangle_vec_within (struct triangle *tri,struct vec *pt);
double triangle_get_y_from_xz(struct triangle *tri,double x, double z);

double triangle_get_min_y(struct triangle* tri);
double triangle_get_max_y(struct triangle* tri);
void ray_tri_flip_y_axis(struct triangle* tri,double y);
void ray_tri_sub_y(struct triangle* tri,double y);
void ray_tri_add_y(struct triangle* tri,double y);

#endif
