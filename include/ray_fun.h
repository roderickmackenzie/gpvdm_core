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

/** @file ray_fun.h
@brief ray tracing header files.
*/
#ifndef ray_fun_h
#define ray_fun_h

#include <vec.h>
#include <sim_struct.h>
#include <device.h>
#include <ray.h>
#define WAIT 0
#define READY 1
#define DONE 2

#define TRUE 1
#define FALSE 0

#define RAY_MAX 5000

void light_update_ray_mat(struct simulation *sim,struct epitaxy *my_epitaxy,struct image *my_image,double lambda);
void image_init(struct image *in);
int between(double v, double x0, double x1);
void add_triangle(struct image *in, double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,int object_uid,int edge);
void ray_reset(struct image *in);
int add_ray(struct simulation *sim,struct image *in,struct vec *start,struct vec *dir,double mag);
void ray_populate_with_shapes(struct simulation *sim,struct image *my_image,struct epitaxy *in);
void dump_plane_to_file(struct simulation *sim,char *file_name,struct image *in);
void dump_plane(struct simulation *sim,struct image *in);
double get_rand();
void obj_norm(struct vec *ret,struct triangle *my_obj);
int ray_intersect(struct vec *ret,struct triangle *my_obj,struct ray *my_ray);
int search_obj(struct simulation *sim,struct image *in,struct ray *my_ray);
int search_triangle(struct simulation *sim,struct image *in,struct ray *my_ray);
int activate_rays(struct image *in);
int pnpoly(struct image *in, struct vec *xy,int id);
void get_refractive(struct simulation *sim,struct image *in,double *alpha,double *n0,double *n1,struct ray *my_ray);
int propergate_next_ray(struct simulation *sim,struct image *in);
struct object *add_box(struct image *in,double x0,double y0,double z0,double dx,double dy,double dz,int sim_edge);
struct object *add_pyramid(struct image *in,double x0,double y0,double z0,double dx,double dy,double dz);
double get_eff(struct image *in);
void light_setup_ray(struct simulation *sim,struct device *cell,struct image *my_image,struct epitaxy *my_epitaxy);
void ray_free(struct simulation *sim,struct image *my_image);
void ray_read_config(struct simulation *sim,struct image *my_image);
void ray_solve(struct simulation *sim,struct device *in, int l);
void ray_solve_all(struct simulation *sim,struct device *in);
void dump_extraction_efficiency(struct simulation *sim,struct image *in);
void dump_ang_escape(struct simulation *sim,struct image *in);
double ray_cal_escape_angle(struct image *in,int l);
void ray_escape_angle_reset(struct image *in,int l);
void ray_load_emission(struct simulation *sim,struct image *my_image);
void dump_ang_escape_as_rgb(struct simulation *sim,struct image *in);
int search_object(struct simulation *sim,struct image *in,struct ray *my_ray);
void ray_malloc(struct simulation *sim,struct image *my_image);
void ray_escape_angle_norm(struct image *in);
double ray_get_avg_extract_eff(struct image *in);
double ray_tri_get_min_y(struct triangle* tri);
double ray_obj_get_min_y(struct simulation *sim,struct image *my_image,struct object *obj);
void ray_object_flip_y_axis(struct simulation *sim,struct image *in,struct object *obj);
void ray_object_sub_y(struct simulation *sim,struct image *in,struct object *obj,double y);
void ray_object_add_y(struct simulation *sim,struct image *in,struct object *obj,double y);

#endif
