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

/** @file triangles.h
@brief Functions to manipulate lots of triangles
*/
#ifndef triangles_io_h
#define triangles_io_h

#include <vec.h>
#include <triangle.h>
void triangle_load_from_file(struct simulation *sim,struct triangles *in,char *file_name);
void triangle_print(struct triangle *in);
void triangles_print(struct triangles *in);
void triangles_free(struct triangles *in);
void triangles_cpy(struct triangles *out,struct triangles *in);
void triangles_find_min(struct vec *out,struct triangles *in);
void triangles_find_max(struct vec *out,struct triangles *in);
void triangles_sub_vec(struct triangles *in,struct vec *v);
void triangles_add_vec(struct triangles *in,struct vec *v);
void triangles_div_vec(struct triangles *in,struct vec *v);
void triangles_mul_vec(struct triangles *in,struct vec *v);
void triangles_cal_edges(struct triangles *in);
void triangles_init(struct triangles *tri);
void triangles_malloc(struct triangles *tri);
void triangles_save(char *file_name,struct triangles *in);
void triangles_add_triangle(struct triangles *obj, double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,int uid,int object_type);
void triangles_set_object_type(struct triangles *in,int object_type);
double triangles_interpolate(struct triangles *in,struct vec *p);
#endif
