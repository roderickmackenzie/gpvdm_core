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

/** @file matrix.c
@brief A struct for the matrix solver
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include "memory.h"
#include "md5.h"
#include "cal_path.h"

void matrix_init(struct matrix *mx)
{
	mx->Ti = NULL;
	mx->Tj = NULL;
	mx->Tx = NULL;
	mx->Txz = NULL;
	mx->b = NULL;
	mx->bz = NULL;
	mx->Tdebug = NULL;
	mx->nz= 0;
	mx->nz_max= 0;
	
	mx->M= 0;

	mx->use_cache=FALSE;
	strcpy(mx->hash,"");
	mx->ittr=0;
	mx->complex_matrix=FALSE;
}

void matrix_cache_reset(struct simulation *sim,struct matrix *mx)
{
	strcpy(mx->hash,"");
	mx->ittr=0;
}

void matrix_dump(struct simulation *sim,struct matrix *mx)
{
int i;
	printf("J:\n");
	for (i=0;i<mx->nz;i++)
	{
		printf_log(sim,"%ld %ld %Le\n",mx->Tj[i],mx->Ti[i],mx->Tx[i]);
	}

	printf("b:");
	
	for (i=0;i<mx->M;i++)
	{
		printf_log(sim,"%Le\n",mx->b[i]);
	}

}

void matrix_dump_b(struct simulation *sim,struct matrix *mx)
{
int i;
	if (mx->complex_matrix==TRUE)
	{
		for (i=0;i<mx->M;i++)
		{
			printf_log(sim,"%Le %Le\n",mx->b[i],mx->bz[i]);
		}
	}else
	{
		for (i=0;i<mx->M;i++)
		{
			printf_log(sim,"%Le\n",mx->b[i]);
		}
	}
}

void matrix_dump_J(struct simulation *sim,struct matrix *mx)
{
int i;
	if (mx->complex_matrix==TRUE)
	{
		for (i=0;i<mx->nz;i++)
		{
			printf_log(sim,"%ld %ld %Le %Le\n",mx->Ti[i],mx->Tj[i],mx->Tx[i],mx->Txz[i]);
		}
	}else
	{

		for (i=0;i<mx->nz;i++)
		{
			printf_log(sim,"%ld %ld %Le\n",mx->Ti[i],mx->Tj[i],mx->Tx[i]);
		}
	}
}


long double matrix_cal_error(struct simulation *sim,struct matrix *mx)
{
int i;
long double sum=0.0;
	for (i=0;i<mx->M;i++)
	{
		sum+=fabsl(mx->b[i]);
	}

return sum;
}

int matrix_solve(struct simulation *sim,struct matrix *mx)
{
char out[100];
struct md5 hash;

	if (mx->use_cache==TRUE)
	{
		if (mx->ittr==0)
		{
			md5_init(&hash);
			md5_update(&hash,(char*)mx->Ti,mx->nz*sizeof(int));
			md5_update(&hash,(char*)mx->Tj,mx->nz*sizeof(int));
			md5_update(&hash,(char*)mx->Tx,mx->nz*sizeof(long double));
			md5_update(&hash,(char*)mx->b,mx->M*sizeof(long double));
			md5_to_str(mx->hash,&hash);

			join_path(2,mx->cache_file_path,get_cache_path(sim),mx->hash);

			if (isfile(mx->cache_file_path)==0)
			{
				return 0;
			}
		}

	}

	if (mx->complex_matrix==FALSE)
	{
		(*sim->dll_matrix_solve)(sim,mx->M,mx->nz,mx->Ti,mx->Tj,mx->Tx,mx->b);
	}else
	{
		(*sim->dll_complex_matrix_solve)(sim,mx->M,mx->nz,mx->Ti,mx->Tj,mx->Tx,mx->Txz,mx->b,mx->bz);
	}

mx->ittr++;
return -1;
}

void matrix_malloc(struct simulation *sim,struct matrix *mx)
{
	mx->Ti=malloc(mx->nz*sizeof(int));
	memset(mx->Ti, 0, mx->nz*sizeof(int));

	mx->Tj=malloc(mx->nz*sizeof(int));
	memset(mx->Tj, 0, mx->nz*sizeof(int));

	mx->Tx=malloc(mx->nz*sizeof(long double));
	memset(mx->Tx, 0, mx->nz*sizeof(long double));

	mx->b=malloc(mx->M*sizeof(long double));
	memset(mx->b, 0, mx->M*sizeof(long double));

	if (mx->complex_matrix==TRUE)
	{
		mx->Txz=malloc(mx->nz*sizeof(long double));
		memset(mx->Txz, 0, mx->nz*sizeof(long double));

		mx->bz=malloc(mx->M*sizeof(long double));
		memset(mx->bz, 0, mx->M*sizeof(long double));
	}
}

void matrix_add_nz_item(struct simulation *sim,struct matrix *mx,int x,int y,long double val)
{
	int i;
	//printf("ns=%d %d\n",mx->nz,mx->nz_max);
	for (i=0;i<mx->nz;i++)
	{
		if ((x==mx->Tj[i])&&(y==mx->Ti[i]))
		{
			//printf("found\n");
			mx->Tx[i]+=val;
			return;
		}
	}

	if (mx->nz>=mx->nz_max)
	{
		mx->nz_max+=1000;
		mx->Ti=realloc(mx->Ti,mx->nz_max*sizeof(int));
		mx->Tj=realloc(mx->Tj,mx->nz_max*sizeof(int));
		mx->Tx=realloc(mx->Tx,mx->nz_max*sizeof(long double));
	}

	mx->Ti[mx->nz]=y;
	mx->Tj[mx->nz]=x;
	mx->Tx[mx->nz]=val;
	mx->nz++;

}
void matrix_realloc(struct simulation *sim,struct matrix *mx)
{
	gdouble *dtemp=NULL;
	int *itemp=NULL;
	itemp = realloc(mx->Ti,mx->nz*sizeof(int));
	if (itemp==NULL)
	{
		ewe(sim,"mx->Ti - memory error\n");
	}else
	{
		mx->Ti=itemp;
		memset(mx->Ti, 0, mx->nz*sizeof(int));
	}

	itemp = realloc(mx->Tj,mx->nz*sizeof(int));
	if (itemp==NULL)
	{
		ewe(sim,"mx->Tj - memory error\n");
	}else
	{
		mx->Tj=itemp;
		memset(mx->Tj, 0, mx->nz*sizeof(int));
	}

	dtemp = realloc(mx->Tx,mx->nz*sizeof(long double));
	if (dtemp==NULL)
	{
		ewe(sim,"mx->Tx - memory error\n");
	}else
	{
		mx->Tx=dtemp;
		memset(mx->Tx, 0, mx->nz*sizeof(long double));
	}

	dtemp = realloc(mx->b,mx->M*sizeof(long double));

	if (dtemp==NULL)
	{
		ewe(sim,"in->b - memory error\n");
	}else
	{
		mx->b=dtemp;
		memset(mx->b, 0, mx->M*sizeof(long double));
	}

	if (mx->complex_matrix==TRUE)
	{
		dtemp = realloc(mx->Txz,mx->nz*sizeof(long double));
		if (dtemp==NULL)
		{
			ewe(sim,"mx->Txz - memory error\n");
		}else
		{
			mx->Txz=dtemp;
			memset(mx->Txz, 0, mx->nz*sizeof(long double));
		}

		dtemp = realloc(mx->bz,mx->M*sizeof(long double));

		if (dtemp==NULL)
		{
			ewe(sim,"in->bz - memory error\n");
		}else
		{
			mx->bz=dtemp;
			memset(mx->bz, 0, mx->M*sizeof(long double));
		}

	}
}

void matrix_zero_b(struct simulation *sim,struct matrix *mx)
{
	memset(mx->b, 0, mx->M*sizeof(long double));

	if (mx->complex_matrix==TRUE)
	{
		memset(mx->bz, 0, mx->M*sizeof(long double));
	}
}

void matrix_save(struct simulation *sim,struct matrix *mx)
{
	char cache_file[PATH_MAX];
	FILE *file;

	join_path(2, cache_file,get_cache_path(sim),mx->hash);

	file=fopen(cache_file,"wb");

	fwrite ((char*)mx->Ti,mx->nz*sizeof(int),1,file);
	fwrite ((char*)mx->Tj,mx->nz*sizeof(int),1,file);
	fwrite ((char*)mx->Tx,mx->nz*sizeof(long double),1,file);
	fwrite ((char*)mx->b,mx->M*sizeof(long double),1,file);

	if (mx->complex_matrix==TRUE)
	{
		fwrite ((char*)mx->Txz,mx->nz*sizeof(long double),1,file);
		fwrite ((char*)mx->bz,mx->M*sizeof(long double),1,file);
	}

	fclose(file);

}

int matrix_load(struct simulation *sim,struct matrix *mx)
{
	char cache_file[PATH_MAX];
	FILE *file;

	join_path(2, cache_file,get_cache_path(sim),mx->hash);

	file=fopen(cache_file,"rb");
	if (file==NULL)
	{
		return -1;
	}

	fread((int*)mx->Ti, mx->nz*sizeof(int), 1, file);
	fread((int*)mx->Tj, mx->nz*sizeof(int), 1, file);
	fread((long double *)mx->Tx, mx->nz*sizeof(long double), 1, file);
	fread((long double *)mx->b, mx->M*sizeof(long double), 1, file);

	if (mx->complex_matrix==TRUE)
	{
		fread((long double *)mx->Txz, mx->nz*sizeof(long double), 1, file);
		fread((long double *)mx->bz, mx->M*sizeof(long double), 1, file);
	}
	fclose(file);

	return 0;

}

void matrix_free(struct simulation *sim,struct matrix *mx)
{
	if (mx->Ti!=NULL)
	{
		free(mx->Ti);
	}

	if (mx->Tj!=NULL)
	{
		free(mx->Tj);
	}

	if (mx->Tx!=NULL)
	{
		free(mx->Tx);
	}

	if (mx->b!=NULL)
	{
		free(mx->b);
	}

	if (mx->Txz!=NULL)
	{
		free(mx->Txz);
	}

	if (mx->bz!=NULL)
	{
		free(mx->bz);
	}

	matrix_init(mx);
}
