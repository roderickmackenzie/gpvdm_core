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

/** @file epitaxy_heat.c
 * @brief Do the maths for the heat parts of the epitaxy
*/

#include <string.h>
#include "epitaxy.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include <cal_path.h>
#include <shape.h>
#include <contacts.h>


long double epitaxy_get_heat_problem_start(struct epitaxy *in)
{
int i=0;
gdouble pos=0.0;
for (i=0;i<in->layers;i++)
{

	//printf("%d\n",in->layer[i].solve_optical_problem);
	//getchar();
	if (in->layer[i].solve_thermal_problem==TRUE)
	{
		return pos;
	}
	pos+=in->layer[i].width;

}

return -1;
}


long double epitaxy_get_heat_problem_stop(struct epitaxy *in)
{
int i=0;
gdouble pos=0.0;
int found=FALSE;
for (i=0;i<in->layers;i++)
{

	if (in->layer[i].solve_thermal_problem==TRUE)
	{
		found=TRUE;
	}

	if ((in->layer[i].solve_thermal_problem==FALSE)&&(found==TRUE))
	{
		return pos;
	}

	pos+=in->layer[i].width;
}

if (found==TRUE)
{
	return pos;
}


return -1;
}

