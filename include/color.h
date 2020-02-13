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

/** @file color.h
	@brief Header file for lib which handles cie_colors
*/

#ifndef cie_color_h
#define cie_color_h
#include <sim_struct.h>
#include <device.h>

void color_cie_init(struct simulation *sim);
void color_cie_load(struct simulation *sim);
void color_cie_free(struct simulation *sim);
void wavelength_to_rgb(int *r,int *g,int *b,double wavelength);
void color_cie_cal_XYZ(struct simulation *sim,long double *X,long double *Y,long double *Z,struct istruct *L_input, int input_in_ev);
void color_XYZ_to_rgb(int *R,int *G, int *B,long double X,long double Y, long double Z);
#endif
