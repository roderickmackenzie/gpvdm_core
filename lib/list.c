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

/** @file inp.c
	@brief Input file interface, files can be in .gpvdm files or stand alone files.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zip.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>

#include "list_struct.h"
#include "util.h"
#include "code_ctrl.h"
#include "const.h"
#include <log.h>
#include <cal_path.h>

void list_free(struct list *in)
{
	int i=0;

	for (i=0;i<in->len;i++)
	{
		free(in->names[i]);
	}

	free(in->names);
}

int list_cmp(struct list *in,char *name)
{
	int i=0;

	for (i=0;i<in->len;i++)
	{
		if (strcmp(name,in->names[i])==0)
		{
			return 0;
		}
	}

return -1;
}

void list_dump(struct list *in)
{
	int i=0;

	for (i=0;i<in->len;i++)
	{
		printf("%s\n",in->names[i]);
	}

}
