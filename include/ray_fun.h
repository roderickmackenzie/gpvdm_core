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

void light_update_ray_mat(struct simulation *sim,struct epitaxy *my_epitaxy,struct image *my_image,long double lambda);
void image_init(struct image *in);
int between(double v, double x0, double x1);
void add_plane(struct image *in, double x0,double y0,double x1,double y1,int id,int edge);
void ray_reset(struct image *in);
void add_ray(struct simulation *sim,struct image *in,struct vec *start,struct vec *dir,double mag);
void dump_plane_to_file(char *file_name,struct image *in);
void dump_plane(struct simulation *sim,struct image *in);
double get_rand();
void obj_norm(struct vec *ret,struct plane *my_obj);
int ray_intersect(struct vec *ret,struct plane *my_obj,struct ray *my_ray);
int search_ojb(struct image *in,int ray);
int activate_rays(struct image *in);
int pnpoly(struct image *in, struct vec *xy,int id);
void get_refractive(struct image *in,double *alpha,double *n0,double *n1,int ray);
int propergate_next_ray(struct simulation *sim,struct image *in);
void add_box(struct image *in,double start_x,double start_y,double x_len,double y_len,int n,int sim_edge);
double get_eff(struct image *in);
void light_setup_ray(struct simulation *sim,struct device *cell,struct image *my_image,struct epitaxy *my_epitaxy);
void ray_free(struct simulation *sim,struct image *my_image);
void ray_read_config(struct simulation *sim,struct image *my_image);
void ray_solve(struct simulation *sim,struct device *in, int l);
void dump_extraction_efficiency(struct simulation *sim,struct image *in);
void dump_ang_escape(struct simulation *sim,struct image *in);
void ray_cal_escape_angle(struct image *in,int l);
#endif
