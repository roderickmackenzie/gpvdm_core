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

#include <stdio.h>
#include <ray.h>
#include <ray_fun.h>
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <device.h>
#include <inp.h>
#include <util.h>
#include <buffer.h>
#include <color.h>

/** @file ray_bat.c
	@brief Materials loader for ray plugin
*/



void ray_load_emission(struct simulation *sim,struct image *my_image)
{
	struct buffer buf;
	long double max;
	long double X;
	long double Y;
	long double Z;
	int R;
	int G;
	int B;

	inter_load(sim,&(my_image->input_spectrum),my_image->input_spectrum_file);
	inter_sort(&(my_image->input_spectrum));

	max=inter_get_max(&(my_image->input_spectrum));
	inter_div_long_double(&(my_image->input_spectrum),max);

	int min_pos=0;
	int max_pos=0;
	for (min_pos=0;min_pos<my_image->input_spectrum.len;min_pos++)
	{
		if (my_image->input_spectrum.data[min_pos]>max*0.05)
		{
			break;
		}
	}

	for (max_pos=my_image->input_spectrum.len-1;max_pos>0;max_pos--)
	{
		if (my_image->input_spectrum.data[max_pos]>max*0.05)
		{
			break;
		}
	}

	//printf("%ld %ld\n",my_image->input_spectrum.len,max_pos);
	//getchar();
	inter_chop(&(my_image->input_spectrum),my_image->input_spectrum.x[min_pos], my_image->input_spectrum.x[max_pos]);

	color_cie_cal_XYZ(sim,&X,&Y,&Z,&(my_image->input_spectrum),FALSE);
	color_XYZ_to_rgb(&R,&G,&B,X,Y,Z);

	buffer_init(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1.0;
	buf.y_mul=1e9;
	strcpy(buf.title,"The input emission spectrum");
	strcpy(buf.type,"linegraph");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Probability");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"a.u.");

	sprintf(buf.rgb,"%.2x%.2x%.2x",R,G,B);

	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=1;
	buf.y=my_image->input_spectrum.len;
	buf.z=1;
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,my_image->input_spectrum.x, my_image->input_spectrum.data, my_image->input_spectrum.len);
	buffer_dump_path(sim,"","emission_input.dat",&buf);
	buffer_free(&buf);

}

