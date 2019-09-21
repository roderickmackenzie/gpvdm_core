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

/** @file hash.c
	@brief Hashing function for fast lookup in tables, but the arrays are now linear so you don't need to hash.
*/

#include <stdio.h>
#include "sim.h"
#include "inp.h"
#include "cal_path.h"
#include "list.h"
#include "md5.h"

void hash_dir(struct simulation *sim,char *out)
{
	struct md5 sum;
	md5_init(&sum);

	int i=0;
	char *buffer;
	unsigned int len;
	long l;
	char newcheck[100];
	struct list files;
	struct inp_file inp;
	inp_listdir(sim, get_input_path(sim),&files);


	for (i=0;i<files.len;i++)
	{
		if (is_numbered_file(files.names[i],"dos")==0)
		{
			//printf("%s\n",files.names[i]);
			inp_read_buffer(sim,&buffer, &l,files.names[i]);
			len=(unsigned int)l;
			
			md5_update(&sum,buffer,len);

			free(buffer);
		}

		if (strcmp(files.names[i],"contacts.inp")==0)
		{
			inp_init(sim,&inp);
			inp_load(sim,&inp,files.names[i]);
			inp_replace(sim,&inp,"#contact_voltage0", "");
			inp_replace(sim,&inp,"#contact_voltage1", "");
			inp_replace(sim,&inp,"#contact_voltage2", "");
			inp_replace(sim,&inp,"#contact_voltage3", "");
			inp_replace(sim,&inp,"#contact_voltage4", "");
			inp_replace(sim,&inp,"#contact_voltage5", "");
			inp_replace(sim,&inp,"#contact_voltage6", "");
			inp_replace(sim,&inp,"#contact_voltage7", "");
			md5_update(&sum,inp.data,inp.fsize);
			inp_free_no_save(sim,&inp);
		}
	}


	list_free(&files);

	md5_to_str(out,&sum);
}

