// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
//
// Copyright (C) 2012-2017 Roderick C. I. MacKenzie info at gpvdm dot com
//
// https://www.gpvdm.com
//
// This program is free software; you can redistribute it and/or modify it
// under the terms and conditions of the GNU Lesser General Public License,
// version 2.1, as published by the Free Software Foundation.
//
// This program is distributed in the hope it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.

/** @file opp.c
@brief RPN opps.
*/


#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#include "util.h"
#include "cal_path.h"
#include "gpvdm_const.h"
#include <rpn.h>
#include <log.h>
#include <sim_struct.h>


void add_opp(struct simulation *sim,struct rpn *in,char *name, int prec, int right_left,void *f)
{
	strcpy(in->opps[in->opp_count].name,name);
	in->opps[in->opp_count].prec=prec;
	in->opps[in->opp_count].right_left=right_left;
	in->opps[in->opp_count].f=f;
	in->opp_count++;
}

int is_opp(struct simulation *sim,struct rpn *in,char *val)
{
	int i;
	for (i=0;i<in->opp_count;i++)
	{
		if (strcmp(val,in->opps[i].name)==0)
		{

			return i;
		}
	}
return -1;
}

char* opp_run(struct simulation *sim,struct rpn *in,char *val,char *out,char* a,char* b)
{
	int i;
	for (i=0;i<in->opp_count;i++)
	{
		if (strcmp(val,in->opps[i].name)==0)
		{
			return (in->opps[i].f)(out,a,b);
		}
	}
return "error";
}

int opp_pr(struct simulation *sim,struct rpn *in,char *val)
{
	return in->opps[is_opp(sim,in,val)].prec;
}

int opp_lr(struct simulation *sim,struct rpn *in,char *val)
{
	return in->opps[is_opp(sim,in,val)].right_left;
}



