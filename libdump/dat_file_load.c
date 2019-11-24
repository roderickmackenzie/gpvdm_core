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

/** @file buffer.c
@brief used to save output files to disk with a nice header, so the user knows what was writtne to them
*/

#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "dat_file.h"
#include "const.h"
#include "code_ctrl.h"
#include "cal_path.h"
#include "dump.h"
#include <log.h>
#include <cache.h>

#include <triangle.h>
#include <triangle_io.h>


void dat_file_load_info(struct simulation *sim,struct dat_file *in,char *file_name)
{
char line[10000];
char *token;
char temp2[200];
int i;
FILE *file;
file=fopen(file_name,"r");

if (file == NULL)
{
	ewe(sim,"Can not load the file %s\n",file_name);
}

do
{
	memset(line,0,10000);
	fgets(line, 10000, file);

	if (line[0]=='#')
	{
		if (strcmp_begin(line,"#x ")==0)
		{
			sscanf(line,"%s %d",temp2,&in->x);
		}else
		if (strcmp_begin(line,"#y ")==0)
		{
			sscanf(line,"%s %d",temp2,&in->y);
		}else
		if (strcmp_begin(line,"#z ")==0)
		{
			sscanf(line,"%s %d",temp2,&in->z);
		}else
		if (strcmp_begin(line,"#type ")==0)
		{
			sscanf(line,"%s %s",temp2,in->type);
		}


	}else
	{
		break;
	}


}while(!feof(file));
fclose(file);

}

