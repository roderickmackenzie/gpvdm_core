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
#include "gpvdm_const.h"
#include "light.h"
#include "device.h"
#include "gpvdm_const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "epitaxy.h"
#include "lang.h"
#include "log.h"
#include <cal_path.h>
#include <light_fun.h>


static int unused __attribute__((unused));

/** @file light_config.c
	@brief Loads the light configuration files.
*/



void light_load_config_file(struct simulation *sim,struct light *li)
{

	char path_temp[PATH_MAX];
	char temp_str[100];

	long double temp=0.0;
	struct inp_file inp;
	struct dim_light *dim=&(li->dim);

	li->disable_transfer_to_electrical_mesh=FALSE;

	join_path(2,path_temp,get_output_path(sim),"optical_output");
	remove_dir(sim,path_temp);

	join_path(2,li->config_file,get_output_path(sim),"light.inp");

	printf_log(sim,"%s: %s\n",_("load"),li->config_file);

	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_input_path(sim),"light.inp");

	inp_check(sim,&inp,1.31);

	inp_search_string(sim,&inp,li->suns_spectrum_file,"#sun");

	inp_search_int(sim,&inp,&li->align_mesh,"#alignmesh");

	inp_search_int(sim,&inp,&dim->ylen,"#meshpoints");

	dim->dy=epitaxy_get_optical_length(li->epi)/((long double)dim->ylen);
	//printf("%Le\n",epitaxy_get_optical_length(li->epi));
	//getchar();
	inp_search_string(sim,&inp,temp_str,"#light_illuminate_from");
	if (strcmp(temp_str,"bottom")==0)
	{
		li->flip_field=TRUE;
	}else
	{
		li->flip_field=FALSE;
	}

	inp_search_int(sim,&inp,&dim->llen,"#lpoints");

	inp_search_gdouble(sim,&inp,&li->lstart,"#lstart");

	inp_search_gdouble(sim,&inp,&li->lstop,"#lstop");

	dim->dl=(li->lstop-li->lstart)/((long double)dim->llen);

	inp_search_gdouble(sim,&inp,&(li->electron_eff),"#electron_eff");
	li->electron_eff=fabs(li->electron_eff);

	inp_search_gdouble(sim,&inp,&(li->hole_eff),"#hole_eff");
	li->hole_eff=fabs(li->hole_eff);

	inp_search_gdouble(sim,&inp,&(temp),"#Psun");
	li->Psun=fabs(temp);

	inp_search_string(sim,&inp,li->mode,"#light_model");

	inp_search_gdouble(sim,&inp,&(li->Dphotoneff),"#Dphotoneff");
	li->Dphotoneff=fabs(li->Dphotoneff);

	inp_search_gdouble(sim,&inp,&(li->ND),"#NDfilter");

	inp_search_string(sim,&inp,li->light_file_generation,"#light_file_generation");

	inp_search_gdouble(sim,&inp,&(li->light_file_generation_shift),"#light_file_generation_shift");

	inp_search_string(sim,&inp,li->light_profile,"#light_profile");

	inp_free(sim,&inp);

}

void light_load_config(struct simulation *sim,struct light *li, struct device *dev)
{
	struct dimensions *dim=&(dev->ns.dim);
	li->epi=&(dev->my_epitaxy);
	li->dim.xlen=dim->xlen;
	li->dim.zlen=dim->zlen;
	li->force_update=FALSE;

	light_load_config_file(sim,li);
	light_load_materials(sim,li,dev);
	light_malloc(sim,li);
	light_build_mesh(sim,li,dim);
	light_build_materials_arrays(sim,li,dev);
}
