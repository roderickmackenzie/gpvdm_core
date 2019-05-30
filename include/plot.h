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

/** @file plot.h
@brief gnuplot interface for interactive plotting.
*/


#ifndef plot_h
#define plot_h
void plot_now_excite(struct simulation *sim);
void plot_open(struct simulation *sim);
void plot_now(struct simulation *sim,struct device *in,char *name);
void plot_close(struct simulation *sim);
void plot_replot(struct simulation *sim);
void set_plot_script_dir(char * in);
#endif
