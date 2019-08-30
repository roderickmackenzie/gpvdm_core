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

/** @file perovskite.c
	@brief A perovskite ion solver.
*/

#include <string.h>
#include <stdlib.h>
#include <dump.h>
#include <dos.h>
#include "sim.h"
#include "solver_interface.h"
#include "dat_file.h"
#include "log.h"
#include <cal_path.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <lang.h>
#include <inp.h>
#include <memory.h>
#include <newton_tricks.h>

long double newton_externalv_simple_perovskite(struct simulation *sim,struct device *in,gdouble V)
{
	long double i0;
	return i0;
}

long double newton_externv_perovskite(struct simulation *sim,struct device *in,gdouble Vtot,int usecap)
{
	int i=0;
	int ii=0;
	long double i0;
	long double i0_last=1000.0;
	long double error=0.0;
	long double first_error=0.0;
	for (i==0;ii<10;ii++)
	{
		for (i=0;i<3;i++)
		{
			i0=newton_externv(sim,in,Vtot,usecap);
			error=fabsl(i0-i0_last);
			//			printf_log(sim,"%s %Le %d\n",_("Perovskite ion"),error,ii);printf_log(sim,"%s %Le %d\n",_("Perovskite ion"),error,ii);

			if (error<1e-6)
			{
				break;
			}

			i0_last=i0;
		}
		//getchar();
	}

	printf_log(sim,"%s %Le\n",_("Electrical+perovskite solver f(x)="),error);
	return i0;
}



