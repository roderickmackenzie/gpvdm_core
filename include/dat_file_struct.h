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

/** @file dat_file_struct.h
	@brief Strcutr to hold .dat files before they are written to disk.
*/

#ifndef dat_file_struct_h
#define dat_file_struct_h
#include "advmath.h"
//#include <zip.h>
#include <triangle.h>

struct dat_file
{
char title[100];
char type[100];
gdouble x_mul;
gdouble y_mul;
gdouble z_mul;
gdouble x_offset;
gdouble y_offset;
gdouble z_offset;
gdouble data_mul;
char x_label[100];
char y_label[100];
char z_label[100];
char data_label[100];
char x_units[100];
char y_units[100];
char z_units[100];
char rgb[100];
char data_units[100];
char section_one[100];
char section_two[100];
int logscale_x;
int logscale_y;
int logscale_z;
int logscale_data;
int write_to_zip;
int norm_x_axis;
int norm_y_axis;
long double data_max;
long double data_min;
int x;
int y;
int z;
gdouble time;
gdouble Vexternal;
char *buf;
struct triangle *data;
int len;
int max_len;
char zip_file_name[400];
//struct zip *zip_file;
};

#endif
