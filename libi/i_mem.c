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




/** @file i_mem.c
	@brief Memory management for i.c
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sim_struct.h>

#include "i.h"
#include "util.h"
#include "cal_path.h"
#include "const.h"
#include <log.h>
#include "inp.h"


static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));


/**Initialize a 1D istruct
@param in input istruct
*/
void inter_init(struct simulation *sim,struct istruct* in)
{
in->len=0;
in->max_len=100;
inter_alloc(in,in->max_len);
}


/**Make a copy of one istruct
@param in input istruct
@param output istruct
@param alloc initialize the memory in the output istruct
*/
void inter_copy(struct istruct* in,struct istruct* orig,int alloc)
{
int i;
in->len=orig->len;

if (alloc==TRUE)
{
inter_alloc(in,orig->len);
}

for  (i=0;i<orig->len;i++)
{
	in->x[i]=orig->x[i];
	in->data[i]=orig->data[i];
}


}

/**Free the structure holding the data
@param in The structure holding the data
*/
void inter_free(struct istruct* in)
{
in->len=0;
in->max_len=0;
free(in->x);
free(in->data);
}


void inter_append(struct istruct* in,long double x,long double y)
{
in->x[in->len]=x;
in->data[in->len]=y;
in->len++;

if ((in->max_len-in->len)<10)
{
	in->max_len+=100;
	inter_realloc(in,in->max_len);
}

}


/**Change the size of an allocated istruct
@param in inout istruct
@param len new length
*/
void inter_realloc(struct istruct* in,int len)
{
in->x=(long double *)realloc (in->x,len*sizeof(long double));
in->data=(long double *)realloc (in->data,len*sizeof(long double));
}

/**Allocate a new 1D istruct
@param in input istruct
@param len new length
*/
void inter_new(struct istruct* in,int len)
{
int i;
strcpy(in->name,"new");

in->len=len;


inter_alloc(in,in->len);

for  (i=0;i<in->len;i++)
{
	in->data[i]=0.0;
}

}

/**Allocate istruct as a 2D array
@param in the array to allocate
@param m number of coloums
@param len length of data to store in the array
*/
void inter_alloc(struct istruct* in,int len)
{
in->x=(long double *)malloc (len*sizeof(long double));
in->data=(long double *)malloc (len*sizeof(long double));
}

void inter_init_mesh(struct istruct* in,int len,long double min,long double max)
{
int i;
in->len=len;
inter_alloc(in,in->len);
memset(in->data, 0, in->len*sizeof(long double));
long double pos=min;
long double dx=(max-min)/((long double)in->len);

for (i=0;i<in->len;i++)
{
	in->x[i]=pos;
	pos+=dx;
}

}
