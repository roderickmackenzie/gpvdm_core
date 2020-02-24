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

/** @file light_materials.c
	@brief This loads in any physical spectra for the light model, not alpha/n data is stored in the epitaxy.
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
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
#include <lang.h>
#include "log.h"
#include <cal_path.h>
#include <dat_file.h>
#include <triangles.h>

static int unused __attribute__((unused));

void light_load_materials(struct simulation *sim,struct light *li, struct device *dev)
{
	printf_log(sim,"%s\n",_("load: materials"));
	struct vec my_vec;
	char file_path[PATH_MAX];

	DIR *theFolder;

	struct inp_file inp;

	theFolder = opendir(get_spectra_path(sim));
	if (theFolder==NULL)
	{
		ewe(sim,_("Optical spectra directory not found\n"));
	}
	closedir (theFolder);
	inp_init(sim,&inp);


	join_path(3,file_path,get_spectra_path(sim),li->suns_spectrum_file,"spectra.inp");

	if (isfile(file_path)!=0)
	{
		ewe(sim,"%s: %s\n",_("File not found"),file_path);
	}

	inter_load(sim,&(li->sun_read),file_path);
	inter_sort(&(li->sun_read));

	inter_mod(&(li->sun_read));

	long double Power=inter_intergrate(&(li->sun_read));
	printf_log(sim,"%s %Le Wm^{-2}\n",_("Power density of the optical spectra:"),Power);

	if (strcmp(li->light_profile,"box")!=0)
	{
		join_path(3,file_path,get_shape_path(sim),li->light_profile,"shape.inp");
		triangle_load_from_file(sim,(&li->light_profile_tri),file_path);

		triangles_find_min(&my_vec,&li->light_profile_tri);
		triangles_sub_vec(&li->light_profile_tri,&my_vec);
		triangles_find_max(&my_vec,&li->light_profile_tri);
		//printf("%Le\n",li->dim.xlen);
		//printf("%Le\n",li->dim.zlen);
		my_vec.x=my_vec.x/dev->xlen;
		my_vec.z=my_vec.z/dev->zlen;
		triangles_div_vec(&li->light_profile_tri,&my_vec);
		triangles_save("test.dat",&li->light_profile_tri);
		//my_vec.x=0.00042;
		//my_vec.z=0.00042;
		//double mul=triangles_interpolate(&li->light_profile_tri,&my_vec);

		//printf("%le %le %le\n",my_vec.x,my_vec.z,mul);
		//getchar();
	}
}


