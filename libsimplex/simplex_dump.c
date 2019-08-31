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

/** @file simplex_dump.c
@brief Dump simplex status
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <simplex.h>
#include <log.h>
#include <sim_struct.h>

void multimin_dump(struct simulation *sim,struct multimin *data)
{
#ifdef simplex_verbose
	int s;
	int d;

	/*	FILE *out;
		out=fopen("one.dat","w");
		fprintf(out,"0.0 0.0\n\n");

		fprintf(out,"50.0 50.0\n\n");

		fprintf(out,"%f %f\n\n",data->center[0],data->center[1]);


		for (s=0;s<data->nsimplex;s++)
		{
			fprintf(out,"%f %f\n",data->p[s][0],data->p[s][1]);
		}
			fprintf(out,"%f %f\n",data->p[0][0],data->p[0][1]);
		fclose(out);

	printf("\n");*/

		for (d=0;d<data->ndim;d++)
		{
			for (s=0;s<data->nsimplex;s++)
			{
				printf_log(sim,"%e ",data->p[s][d]);
			}

			printf_log(sim,"\n");
		}

		for (s=0;s<data->nsimplex;s++)
		{
			printf_log(sim,"y[%d] = %e \n",s,data->y[s]);
		}

		for (d=0;d<data->ndim;d++)
		{
			printf_log(sim,"center[%d] = %e \n",d,data->center[d]);
		}


		printf_log(sim,"0worst=%d 1worst=%d best='%d'\n",data->i_hi0,data->i_hi1,data->i_lo);

		printf_log(sim,"error=%e\n",data->error);

	//printf("\n");

	//getchar();
#endif
}
