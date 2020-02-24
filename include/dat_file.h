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

/** @file dat_file.h
	@brief Strcutr to hold .dat files before they are written to disk.
*/

#ifndef dat_file_h
#define dat_file_h
#include "advmath.h"
#include <dat_file_struct.h>
//#include <zip.h>
#include <device.h>
#include <triangle.h>


void buffer_zip_set_name(struct dat_file *in,char * name);
void buffer_init(struct dat_file *in);
void buffer_malloc(struct dat_file *in);
void buffer_add_xy_data(struct simulation *sim,struct dat_file *in,gdouble *x, gdouble *y, int len);
void buffer_add_string(struct dat_file *in,char * string);
void buffer_add_info(struct simulation *sim,struct dat_file *in);
void buffer_dump(struct simulation *sim,char * file,struct dat_file *in);
void buffer_dump_path(struct simulation *sim,char *path,char * file,struct dat_file *in);
void buffer_free(struct dat_file *in);
void buffer_dump_aes(char *path,char * file,struct dat_file *in,char *key_text);
void buffer_add_xy_data_z_label(struct dat_file *in,gdouble *x, gdouble *y, gdouble *z, int len);
void buffer_add_zxy_long_double_light_data(struct simulation *sim,struct dat_file *in,long double ***data, struct dim_light *dim);
void buffer_add_zxy_heat_data(struct simulation *sim,struct dat_file *in,long double ***data, struct dim_heat *dim);
void buffer_dump_cache(struct simulation *sim,char * file,struct dat_file *in);
void buffer_add_dir(struct simulation *sim,char * file_name);
void buffer_add_3d_device_data_including_boundaries(struct simulation *sim,struct dat_file *buf,struct device *in,gdouble ***data,long double **left,long double **right);
void buffer_add_2d_device_data_int(struct simulation *sim,struct dat_file *buf,struct device *in,int **data);

void dat_file_load_info(struct simulation *sim,struct dat_file *in,char *file_name);
#endif
