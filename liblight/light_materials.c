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
#include <lang.h>
#include "log.h"
#include <cal_path.h>
#include <dat_file.h>

static int unused __attribute__((unused));

void light_load_materials(struct simulation *sim,struct light *in)
{
printf_log(sim,"%s\n",_("load: materials"));
char file_path[PATH_MAX];

DIR *theFolder;

struct inp_file inp;

/////////////////////////////////////////////////////
theFolder = opendir(get_spectra_path(sim));
if (theFolder==NULL)
{
	ewe(sim,_("Optical spectra directory not found\n"));
}
closedir (theFolder);
inp_init(sim,&inp);


	join_path(3,file_path,get_spectra_path(sim),in->suns_spectrum_file,"spectra.inp");

	if (isfile(file_path)!=0)
	{
		ewe(sim,"%s: %s\n",_("File not found"),file_path);
	}

inter_load(sim,&(in->sun_read),file_path);
inter_sort(&(in->sun_read));

inter_mod(&(in->sun_read));

long double Power=inter_intergrate(&(in->sun_read));
printf_log(sim,"%s %Le Wm^{-2}\n",_("Power density of the optical spectra:"),Power);
}


