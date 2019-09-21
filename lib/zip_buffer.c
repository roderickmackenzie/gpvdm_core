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

/** @file zip_buffer.c
	@brief Writing and reading zip buffer
*/

#include <stdio.h>
#include <stdlib.h>

	#include <zlib.h>

#include <code_ctrl.h>
#include <sim.h>
#include <inp.h>
#include <util.h>
#include <dat_file.h>
#include <epitaxy.h>
#include <lang.h>

static int unused __attribute__((unused));

#include "dump.h"
#include "log.h"
#include "cal_path.h"


void write_zip_buffer(struct simulation *sim,char *outfile,long double *buf,int buf_len)
{
		gzFile file;
		file = gzopen (outfile, "w9b");
		gzwrite (file, (char*)buf, buf_len*sizeof(gdouble));
		gzclose (file);
	FILE * yes;
	yes = fopen (outfile, "ab");
	int temp1=buf_len*sizeof(gdouble);
	fwrite ((char*)&temp1, sizeof(int),1,yes);
	fclose (yes);

	//out = fopen(outfile, "wb");
	//fwrite((char*)buf, buf_len*sizeof(gdouble), 1, out);
	//fclose(out);

}


int read_zip_buffer(struct simulation *sim,char *file_name,long double **buf)
{
	int len;

	FILE *tl=fopen(file_name,"rb");
	if (tl==NULL)
	{
		return -1;
	}

	fseek(tl, -4, SEEK_END);
	if (fread((char*)&len, 4, 1, tl)==0)
	{
		return -1;
	}

	fclose(tl);

		gzFile file_in;
		file_in = gzopen (file_name, "rb");
		if (file_in==Z_NULL)
		{
			ewe(sim,_("File not found\n"));
		}



	int buf_len=len/sizeof(long double);

	(*buf)=(long double *)malloc(sizeof(long double)*buf_len);


		gzread (file_in, (char*)(*buf), len);
		gzclose(file_in);

	return buf_len;
}


