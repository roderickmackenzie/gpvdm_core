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

/** @file epitaxy_materials.c
	@brief Load the materials into the epitaxy structure.
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include "util.h"
#include "const.h"
#include "device.h"
#include "const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "epitaxy.h"
#include <lang.h>
#include "log.h"
#include <cal_path.h>
#include <dat_file.h>
#include <color.h>



static int unused __attribute__((unused));

void epitaxy_load_emission(struct simulation *sim,struct epi_layer *layer)
{
	char temp_file[100];
	struct dat_file buf;
	long double max;
	long double X;
	long double Y;
	long double Z;
	int R;
	int G;
	int B;
	int max_pos;

	if (inp_isfile(sim,layer->pl_spectrum_file)!=0)
	{
		ewe(sim,"The emission file %s does not exist",layer->pl_spectrum_file);
	}

	inter_load(sim,&(layer->pl_spectrum),layer->pl_spectrum_file);
	inter_sort(&(layer->pl_spectrum));

	max=inter_get_max(&(layer->pl_spectrum));

	inter_div_long_double(&(layer->pl_spectrum),max);
	inter_mul(&(layer->pl_spectrum),layer->pl_experimental_emission_efficiency);

	max_pos=inter_get_max_pos(&(layer->pl_spectrum));
	layer->peak_wavelength=layer->pl_spectrum.x[max_pos];


	int min_pos=0;
	for (min_pos=0;min_pos<layer->pl_spectrum.len;min_pos++)
	{
		if (layer->pl_spectrum.data[min_pos]>max*0.05)
		{
			break;
		}
	}


	for (max_pos=layer->pl_spectrum.len-1;max_pos>0;max_pos--)
	{
		if (layer->pl_spectrum.data[max_pos]>max*0.05)
		{
			break;
		}
	}

	//printf("%ld %ld\n",layer->pl_spectrum.len,max_pos);
	//getchar();
	inter_chop(&(layer->pl_spectrum),layer->pl_spectrum.x[min_pos], layer->pl_spectrum.x[max_pos]);

	color_cie_cal_XYZ(sim,&X,&Y,&Z,&(layer->pl_spectrum),FALSE);

	color_XYZ_to_rgb(&R,&G,&B,X,Y,Z);

	buffer_init(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1.0;
	buf.y_mul=1e9;

	sprintf(buf.title,"The input emission spectrum for layer %s",layer->name);

	strcpy(buf.type,"linegraph");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Probability");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"a.u.");


	sprintf(buf.rgb,"%.2x%.2x%.2x",R,G,B);

	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=1;
	buf.y=layer->pl_spectrum.len;
	buf.z=1;
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,layer->pl_spectrum.x, layer->pl_spectrum.data, layer->pl_spectrum.len);
	sprintf(temp_file,"emission_input_%d.dat",layer->layer_number);
	buffer_dump_path(sim,"",temp_file,&buf);
	buffer_free(&buf);
	printf("e\n");

}

void epitaxy_load_materials(struct simulation *sim,struct epitaxy *in)
{
	printf_log(sim,"%s\n",_("epitaxy load: materials"));
	struct dat_file buf;
	int i=0;
	char temp[100];
	char fit_file[1000];
	char file_path[400];

	DIR *theFolder;

	struct inp_file inp;

	theFolder = opendir(get_materials_path(sim));
	if (theFolder==NULL)
	{
		ewe(sim,_("No optical materials found\n"));
	}
	closedir (theFolder);


	gdouble alpha_mul=1.0;
	gdouble n_mul=1.0;
	gdouble wavelength_shift_n=0.0;
	gdouble wavelength_shift_alpha=0.0;

	inp_init(sim,&inp);
	char patch_file[400];
	char out_file[400];
	char token[400];
	int ii=0;
	gdouble b=0.0;
	gdouble a=0.0;
	gdouble c=0.0;
	char type[40];
	//printf("load materails\n");
	//getchar();

	for (i=0;i<in->layers;i++)
	{

		//mat.inp
		join_path(3, fit_file,get_materials_path(sim),in->mat_file[i],"mat.inp");

		inp_load(sim,&inp,fit_file);

		char default_file_alpha[100];
		char default_file_n[100];

		inp_search_string(sim,&inp,default_file_alpha,"#mat_default_file_alpha");
		inp_search_string(sim,&inp,default_file_n,"#mat_default_file_n");
		inp_get_array_gdouble(sim,in->rgb[i],&inp,"#red_green_blue");
		inp_free(sim,&inp);

		join_path(3, file_path,get_materials_path(sim),in->mat_file[i],"alpha.gmat");

		if (isfile(file_path)!=0)
		{
			ewe(sim,"%s: %s\n",_("File not found"),file_path);
		}

		inter_load(sim,&(in->layer[i].alpha),file_path);
		inter_sort(&(in->layer[i].alpha));

		join_path(3,file_path,get_materials_path(sim),in->mat_file[i],"n.gmat");

		if (isfile(file_path)!=0)
		{
			ewe(sim,"%s: %s\n",_("File not found"),file_path);
		}

		inter_load(sim,&(in->layer[i].n),file_path);
		inter_sort(&(in->layer[i].n));

		//printf("n %d\n",in->layer[i].n.len);
		//printf("alpha %d\n",in->layer[i].alpha.len);
		//getchar();
		//inter_dump(sim,&in->layer[i].n);
		//printf_log(sim,"%s %d\n",file_path,in->layer[i].n.len);

		//getchar();

	}

}


void epitaxy_free_materials(struct epitaxy *in)
{
int i;
	for (i=0;i<in->layers;i++)
	{
		inter_free(&(in->layer[i].n));
		inter_free(&(in->layer[i].alpha));
		if (strcmp(in->layer[i].pl_spectrum_file,"none")!=0)
		{
			inter_free(&(in->layer[i].pl_spectrum));
		}
	}

}

