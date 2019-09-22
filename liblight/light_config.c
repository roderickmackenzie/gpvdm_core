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



#include "util.h"
#include "const.h"
#include "light.h"
#include "device.h"
#include "const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "epitaxy.h"
#include "lang.h"
#include "log.h"
#include <cal_path.h>

static int unused __attribute__((unused));

/** @file light_config.c
	@brief Loads the light configuration files.
*/

void light_free_epitaxy(struct light *in)
{
	free(in->G_percent);
}

void light_import_epitaxy(struct simulation *sim,struct light *in,struct epitaxy *my_epitaxy)
{
	int i;
	in->my_epitaxy=my_epitaxy;
	in->force_update=FALSE;
	

	in->G_percent=(gdouble *)malloc(my_epitaxy->layers*sizeof(gdouble));

	in->ylen=0.0;
	in->device_ylen=0.0;

	for (i=0;i<my_epitaxy->layers;i++)
	{
		in->ylen+=my_epitaxy->layer[i].width;
	}

in->device_start=epitaxy_get_device_start(my_epitaxy);

if (in->device_start!=-1)
{
	in->device_start_layer=epitaxy_get_device_start_i(my_epitaxy);
	in->device_ylen=epitaxy_get_electrical_length(my_epitaxy);
}else
{
	in->device_start=0.0;
	in->device_start_layer=0;
	in->device_ylen=0.0;
}

}

void light_load_config_file(struct simulation *sim,struct light *in)
{
	char path_temp[1000];
	char temp_str[100];

	gdouble temp=0.0;
	struct inp_file inp;

	in->disable_transfer_to_electrical_mesh=FALSE;

	join_path(2,path_temp,get_output_path(sim),"optical_output");
	remove_dir(sim,path_temp);

	join_path(2,in->config_file,get_output_path(sim),"light.inp");

	printf_log(sim,"%s: %s\n",_("load"),in->config_file);

	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_input_path(sim),"light.inp");

	inp_check(sim,&inp,1.31);

	inp_search_string(sim,&inp,in->suns_spectrum_file,"#sun");

	inp_search_int(sim,&inp,&in->align_mesh,"#alignmesh");

	inp_search_int(sim,&inp,&in->points,"#meshpoints");

	inp_search_string(sim,&inp,temp_str,"#light_illuminate_from");
	if (strcmp(temp_str,"bottom")==0)
	{
		in->flip_field=TRUE;
	}else
	{
		in->flip_field=FALSE;
	}

	inp_search_int(sim,&inp,&in->lpoints,"#lpoints");

	inp_search_gdouble(sim,&inp,&in->lstart,"#lstart");

	inp_search_gdouble(sim,&inp,&in->lstop,"#lstop");

	inp_search_gdouble(sim,&inp,&(in->Eg),"#Eg");

	inp_search_gdouble(sim,&inp,&(in->electron_eff),"#electron_eff");
	in->electron_eff=fabs(in->electron_eff);

	inp_search_gdouble(sim,&inp,&(in->hole_eff),"#hole_eff");
	in->hole_eff=fabs(in->hole_eff);

	inp_search_gdouble(sim,&inp,&(temp),"#Psun");
	in->Psun=fabs(temp);

	inp_search_string(sim,&inp,in->mode,"#light_model");

	inp_search_gdouble(sim,&inp,&(in->Dphotoneff),"#Dphotoneff");
	in->Dphotoneff=fabs(in->Dphotoneff);

	inp_search_gdouble(sim,&inp,&(in->ND),"#NDfilter");

	inp_search_gdouble(sim,&inp,&(temp),"#high_sun_scale");

	inp_search_string(sim,&inp,in->light_file_generation,"#light_file_generation");

	in->Psun*=fabs(temp);

	inp_search_gdouble(sim,&inp,&(in->light_file_generation_shift),"#light_file_generation_shift");

	inp_free(sim,&inp);

}

