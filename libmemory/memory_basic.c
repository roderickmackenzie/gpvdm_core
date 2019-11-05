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

/** @file memory_basic.c
@brief memory functions for 3D arrays
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include <complex_solver.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include "memory.h"

void free_srh_bands(struct dimensions *dim, gdouble **** var)
{
	if (var==NULL)
	{
		return;
	}

	int x=0;
	int y=0;
	int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
					free(var[z][x][y]);
			}
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

	*var=NULL;
}

void malloc_3d_gdouble(struct dimensions *dim, gdouble * (***var))
{
	int x=0;
	int y=0;
	int z=0;


	*var = (gdouble ***) malloc(dim->zmeshpoints * sizeof(gdouble **));

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		(*var)[z] = (gdouble **) malloc(dim->xmeshpoints * sizeof(gdouble*));
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			(*var)[z][x] = (gdouble *) malloc(dim->ymeshpoints * sizeof(gdouble));
			memset((*var)[z][x], 0, dim->ymeshpoints * sizeof(gdouble));
		}
	}

}

void three_d_set_gdouble(struct dimensions *dim, gdouble ***var, gdouble val)
{
int x=0;
int y=0;
int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				var[z][x][y]=val;
			}

		}
	}

}


void three_d_sub_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***sub)
{
int x=0;
int y=0;
int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				var[z][x][y]-=sub[z][x][y];
			}

		}
	}

}

void three_d_copy_gdouble(struct dimensions *dim, gdouble ***dst, gdouble ***src)
{
int x=0;
int y=0;
int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				dst[z][x][y]=src[z][x][y];
			}

		}
	}

}

void srh_copy_gdouble(struct dimensions *dim, gdouble ****dst, gdouble ****src)
{
int x=0;
int y=0;
int z=0;
int b=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				for (b = 0; b < dim->srh_bands; b++)
				{
					dst[z][x][y][b]=src[z][x][y][b];
				}
			}

		}
	}

}

void three_d_add_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***add)
{
int x=0;
int y=0;
int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				var[z][x][y]+=add[z][x][y];
			}

		}
	}

}

void three_d_mul_gdouble(struct dimensions *dim, gdouble ***src, gdouble val)
{
int x=0;
int y=0;
int z=0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				src[z][x][y]*=val;
			}

		}
	}

}

long double three_d_avg(struct device *in, long double ***src)
{
int x=0;
int y=0;
int z=0;
long double sum=0.0;
long double ret=0.0;

long double dx=0.0;
long double dy=0.0;
long double dz=0.0;

struct dimensions *dim=&(in->ns.dim);

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				if (y==0)

				sum+=src[z][x][y]*dim->dxmesh[x]*dim->dymesh[y]*dim->dzmesh[z];
//				printf("%Le %Le %Le %Le %Le %Le\n",dim->dxmesh[x],dim->dymesh[y],dim->dzmesh[z],in->zlen,in->xlen,in->ylen);
			}

		}
	}

ret=sum/(in->zlen*in->xlen*in->ylen);
//printf("ret=%Le\n",ret);
return ret;
}

void three_d_printf(struct dimensions *dim, long double ***src)
{
int x=0;
int y=0;
int z=0;
long double sum=0.0;
long double ret=0.0;
	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				printf("%Le\n",src[z][x][y]);
			}

		}
	}

return;
}

long double three_d_avg_fabsl(struct device *in, long double ***src)
{
int x=0;
int y=0;
int z=0;
long double sum=0.0;
long double ret=0.0;
struct dimensions *dim=&(in->ns.dim);
	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				sum+=fabsl(src[z][x][y])*dim->dxmesh[x]*dim->dymesh[y]*dim->dzmesh[z];
			}

		}
	}

ret=sum/(in->zlen*in->xlen*in->ylen);
return ret;
}

long double three_d_integrate(struct dimensions *dim, long double ***src)
{
int x=0;
int y=0;
int z=0;
long double sum=0.0;

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				sum+=src[z][x][y]*dim->dxmesh[x]*dim->dymesh[y]*dim->dzmesh[z];
			}

		}
	}

return sum;
}

void free_3d_gdouble(struct dimensions *dim, gdouble ***var)
{
	int x=0;
	int y=0;
	int z=0;

	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zmeshpoints; z++)
	{

		for (x = 0; x < dim->xmeshpoints; x++)
		{
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

	*var=NULL;

}

void malloc_zx_gdouble(struct dimensions *dim, gdouble * (**var))
{
	int z=0;

	*var = (gdouble **) malloc(dim->zmeshpoints * sizeof(gdouble *));

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		(*var)[z] = (gdouble *) malloc(dim->xmeshpoints * sizeof(gdouble));
		memset((*var)[z], 0, dim->xmeshpoints * sizeof(gdouble));
	}

}

void mem_set_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in)
{
	int z=0;
	int x=0;


	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			data_out[z][x]=data_in[z][x];
		}
	}

}

void mem_add_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in)
{
	int z=0;
	int x=0;


	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			data_out[z][x]+=data_in[z][x];
		}
	}

}
void free_zx_gdouble(struct dimensions *dim, gdouble **var)
{
	int z=0;

	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		free(var[z]);
	}

	free(var);

	*var=NULL;
}

void malloc_zx_int(struct dimensions *dim, int * (**var))
{
	int z=0;

	*var = (int **) malloc(dim->zmeshpoints * sizeof(int *));

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		(*var)[z] = (int *) malloc(dim->xmeshpoints * sizeof(int));
		memset((*var)[z], 0, dim->xmeshpoints * sizeof(int));
	}

}

void free_zx_int(struct dimensions *dim, int **var)
{
	int z=0;

	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		free(var[z]);
	}

	free(var);

	*var=NULL;
}


void malloc_3d_int(struct dimensions *dim, int * (***var))
{
	int x=0;
	int y=0;
	int z=0;


	*var = (int ***) malloc(dim->zmeshpoints * sizeof(int **));

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		(*var)[z] = (int **) malloc(dim->xmeshpoints * sizeof(int*));
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			(*var)[z][x] = (int *) malloc(dim->ymeshpoints * sizeof(int));
			memset((*var)[z][x], 0, dim->ymeshpoints * sizeof(int));
		}
	}

}

void free_3d_int(struct dimensions *dim, int ***var)
{
	int x=0;
	int y=0;
	int z=0;

	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zmeshpoints; z++)
	{

		for (x = 0; x < dim->xmeshpoints; x++)
		{
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

}

void memory_flip_1d_int(int *var,int len)
{
	int x=0;
	int y=0;
	int z=0;
	int * data=malloc(sizeof(int)*len);
	for (y=0;y<len;y++)
	{
		data[y]=var[len-1-y];
	}

	for (y=0;y<len;y++)
	{
		var[y]=data[y];
	}

	free(data);

}

void memory_flip_1d_long_double(long double *var,int len)
{
	int x=0;
	int y=0;
	int z=0;
	long double * data=malloc(sizeof(long double)*len);
	for (y=0;y<len;y++)
	{
		data[y]=var[len-1-y];
	}

	for (y=0;y<len;y++)
	{
		var[y]=data[y];
	}

	free(data);

}
void malloc_srh_bands(struct dimensions *dim, gdouble * (****var))
{
	int x=0;
	int y=0;
	int z=0;

	//printf("alloc %d %d %d %d \n",dim->xmeshpoints,dim->ymeshpoints,dim->zmeshpoints,dim->srh_bands);

	*var = (gdouble ****) malloc(dim->zmeshpoints * sizeof(gdouble ***));

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		(*var)[z] = (gdouble ***) malloc(dim->xmeshpoints * sizeof(gdouble**));
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			(*var)[z][x] = (gdouble **) malloc(dim->ymeshpoints * sizeof(gdouble*));
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				if (dim->srh_bands != 0)
				{
					(*var)[z][x][y] = (gdouble *) malloc(dim->srh_bands * sizeof(gdouble));
					memset((*var)[z][x][y], 0, dim->srh_bands * sizeof(gdouble));
				}else
				{
					(*var)[z][x][y] = NULL;
				}
			}

		}
	}

}

void three_d_interpolate_gdouble(long double ***out, long double ***in, struct dimensions *dim_out, struct dimensions *dim_in)
{
int x=0;
int y=0;
int z=0;

int yi;
int xi;

long double y_out;
long double x_out;

long double y00;
long double y01;
long double yr;
long double y0;

long double y10;
long double y11;
long double y1;

long double x0;
long double x1;
long double xr;

long double c;

	z=0;
	for (x = 0; x < dim_out->xmeshpoints; x++)
	{

		x_out=dim_out->xmesh[x];
		xi=hashget(dim_in->xmesh,dim_in->xmeshpoints,x_out);

		for (y = 0; y < dim_out->ymeshpoints; y++)
		{
			y_out=dim_out->ymesh[y];
			yi=hashget(dim_in->ymesh,dim_in->ymeshpoints,y_out);

			y00=dim_in->ymesh[yi];
			y01=dim_in->ymesh[yi+1];
			yr=(y_out-y00)/(y01-y00);
			y0=in[z][xi][yi]+yr*(in[z][xi][yi+1]-in[z][xi][yi]);

			y10=dim_in->ymesh[yi];
			y11=dim_in->ymesh[yi+1];
			yr=(y_out-y10)/(y11-y10);
			y1=in[z][xi+1][yi]+yr*(in[z][xi+1][yi+1]-in[z][xi+1][yi]);

			x0=dim_in->xmesh[xi];
			x1=dim_in->xmesh[xi+1];
			xr=(x_out-x0)/(x1-x0);

			c=y0+xr*(y1-y0);
			out[z][x][y]=c;
		}

	}

}

void three_d_interpolate_srh(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band)
{
int x=0;
int y=0;
int z=0;

int yi;
int xi;

long double y_out;
long double x_out;

long double y00;
long double y01;
long double yr;
long double y0;

long double y10;
long double y11;
long double y1;

long double x0;
long double x1;
long double xr;

long double c;

	z=0;
	for (x = 0; x < dim_out->xmeshpoints; x++)
	{

		x_out=dim_out->xmesh[x];
		xi=hashget(dim_in->xmesh,dim_in->xmeshpoints,x_out);

		for (y = 0; y < dim_out->ymeshpoints; y++)
		{
			y_out=dim_out->ymesh[y];
			yi=hashget(dim_in->ymesh,dim_in->ymeshpoints,y_out);

			y00=dim_in->ymesh[yi];
			y01=dim_in->ymesh[yi+1];
			yr=(y_out-y00)/(y01-y00);
			y0=in[z][xi][yi][band]+yr*(in[z][xi][yi+1][band]-in[z][xi][yi][band]);

			y10=dim_in->ymesh[yi];
			y11=dim_in->ymesh[yi+1];
			yr=(y_out-y10)/(y11-y10);
			y1=in[z][xi+1][yi][band]+yr*(in[z][xi+1][yi+1][band]-in[z][xi+1][yi][band]);

			x0=dim_in->xmesh[xi];
			x1=dim_in->xmesh[xi+1];
			xr=(x_out-x0)/(x1-x0);

			c=y0+xr*(y1-y0);
			out[z][x][y][band]=c;
		}

	}

}

void three_d_interpolate_srh2(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band)
{
int x=0;
int y=0;
int z=0;

int yi;
int xi;

long double y_out;
long double x_out;

long double y00;
long double y01;
long double yr;
long double y0;

long double y10;
long double y11;
long double y1;

long double x0;
long double x1;
long double xr;

long double c;

	z=0;
	for (x = 0; x < dim_out->xmeshpoints; x++)
	{

		x_out=dim_out->xmesh[x];
		xi=hashget(dim_in->xmesh,dim_in->xmeshpoints,x_out);

		for (y = 0; y < dim_out->ymeshpoints; y++)
		{
			y_out=dim_out->ymesh[y];
			yi=hashget(dim_in->ymesh,dim_in->ymeshpoints,y_out);

			y00=dim_in->ymesh[yi];
			y01=dim_in->ymesh[yi+1];
			yr=(y_out-y00)/(y01-y00);
			y0=in[z][xi][yi][band]+yr*(in[z][xi][yi+1][band]-in[z][xi][yi][band]);

			y10=dim_in->ymesh[yi];
			y11=dim_in->ymesh[yi+1];
			yr=(y_out-y10)/(y11-y10);
			y1=in[z][xi+1][yi][band]+yr*(in[z][xi+1][yi+1][band]-in[z][xi+1][yi][band]);

			x0=dim_in->xmesh[xi];
			x1=dim_in->xmesh[xi+1];
			xr=(x_out-x0)/(x1-x0);

			c=y0+xr*(y1-y0);
			out[z][x][y][band]=c;
		}

	}

}

void three_d_quick_dump(char *file_name, long double ***in, struct dimensions *dim)
{
int x=0;
int y=0;
int z=0;
	FILE *out=fopen(file_name,"w");

	for (z = 0; z < dim->zmeshpoints; z++)
	{

		for (x = 0; x < dim->xmeshpoints; x++)
		{

			for (y = 0; y < dim->ymeshpoints; y++)
			{
				fprintf(out,"%Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],in[z][x][y]);
			}

			fprintf(out,"\n");
		}
	}

fclose(out);
}

void srh_quick_dump(char *file_name, long double ****in, struct dimensions *dim,int band)
{
int x=0;
int y=0;
int z=0;
	FILE *out=fopen(file_name,"w");

	for (z = 0; z < dim->zmeshpoints; z++)
	{

		for (x = 0; x < dim->xmeshpoints; x++)
		{

			for (y = 0; y < dim->ymeshpoints; y++)
			{
				fprintf(out,"%Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],in[z][x][y][band]);
			}

			fprintf(out,"\n");
		}
	}

fclose(out);
}
