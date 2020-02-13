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

/** @file checksum.c
	@brief An MD5 like check sum for deciding if files have been changed.  Implemented because I do not want openssl dependencies.
*/

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "checksum.h"
#include <gpvdm_const.h>
#include "inp.h"
#include "util.h"
#include "cal_path.h"
#include "md5.h"
static int unused __attribute__((unused));

void checksum_write(struct simulation *sim,char *file_name)
{

FILE *file;
char *buffer;
unsigned long len;
long l;
char chkfile[100];
char temp[100];
struct md5 sum;
md5_init(&sum);

get_file_name_from_path(temp,file_name);

sprintf(chkfile,"%s.chk",temp);

inp_read_buffer(sim,&buffer, &l,file_name);
len=(unsigned int)l;

md5_update(&sum,buffer,len);
md5_to_str(temp,&sum);

free(buffer);

file=fopen(chkfile,"w");
if (file==NULL)
{
	ewe(sim,"File %s not found\n",chkfile);
}
fprintf(file,"%s\n",temp);
fclose(file);
}

int checksum_check(struct simulation *sim,char *file_name)
{
struct md5 sum;
md5_init(&sum);

FILE *file;
char *buffer;
unsigned long len;
char chkfile[100];
char newcheck[100];
char fromfile[100];

strcpy(newcheck,"hello");

sprintf(chkfile,"%s.chk",file_name);
long l;
inp_read_buffer(sim,&buffer, &l,file_name);
len=(unsigned int)l;

md5_update(&sum,buffer,len);
md5_to_str(newcheck,&sum);
//checksum(newcheck,buffer, len);
free(buffer);

file=fopen(chkfile,"r");

if (!file)
{
	return FALSE;
}

unused=fscanf(file,"%s\n",fromfile);
fclose(file);

if (strcmp(newcheck,fromfile)==0)
{
	return TRUE;
}else
{
	return FALSE;
}

return 0;
}



