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



void zxy_malloc_gdouble(struct dimensions *dim, gdouble * (***var))
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

long double zxy_min_gdouble(struct dimensions *dim, gdouble ***var)
{
	int x=0;
	int y=0;
	int z=0;

	long double min=var[0][0][0];

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				if (var[z][x][y]<min)
				{
					min=var[z][x][y];
				}
			}
		}
	}

return min;
}

long double zxy_max_gdouble(struct dimensions *dim, gdouble ***var)
{
	int x=0;
	int y=0;
	int z=0;

	long double max=var[0][0][0];

	for (z = 0; z < dim->zmeshpoints; z++)
	{
		for (x = 0; x < dim->xmeshpoints; x++)
		{
			for (y = 0; y < dim->ymeshpoints; y++)
			{
				if (var[z][x][y]>max)
				{
					max=var[z][x][y];
				}
			}
		}
	}

return max;
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

void free_3d_gdouble(struct dimensions *dim, gdouble * (***in_var))
{
	int x=0;
	int y=0;
	int z=0;

	long double ***var=*in_var;
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

	*in_var=NULL;

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


void zxy_load_long_double(struct simulation *sim, struct dimensions *dim,long double **** data,char *file_name)
{
	char line[1000];
	FILE *file;
	int data_found=FALSE;
	int items_per_line=0;
	int x;
	int y;
	int z;
	long double x_val;
	long double y_val;
	long double z_val;
	long double val;
	long double ***dat=*data;
	struct dat_file d;
	dat_file_load_info(sim,&d,file_name);
	if ((d.x!=dim->xmeshpoints)||(d.y!=dim->ymeshpoints)||(d.z!=dim->zmeshpoints))
	{
		ewe(sim,"not matching dim\n");
	}

	//zxy_malloc_gdouble(dim, data);

	items_per_line++;

	if (d.y>1)
	{
		items_per_line++;	
	}

	if (d.x>1)
	{
		items_per_line++;	
	}

	if (d.z>1)
	{
		items_per_line++;	
	}

	file=fopen(file_name,"r");
	x=0;
	y=0;
	z=0;
	int ret=0;
	do
	{
		memset(line,0,1000);
		ret=gpvdm_fgets(line, 1000, file);

		if (strcmp(line,"#end")==0)
		{
			break;
		}

		if (data_found==TRUE)
		{
			if (ret>0)
			{

				if (items_per_line==3)
				{
					sscanf(line,"%Le %Le %Le",&x_val,&z_val,&val);
					dat[z][x][y]=val;
					//printf("%Le %Le %Le %d %d %d\n",x_val,z_val,dat[z][x][y],z,x,y);
					//getchar();
				}else
				{
					ewe(sim,"I don't know how to read this type of file\n");
				}

				y++;
				if (y>=dim->ymeshpoints)
				{
					y=0;
					x++;
					if (x>=dim->xmeshpoints)
					{
						x=0;
						z++;
						if (z>=dim->zmeshpoints)
						{
							z=0;
						}
					}
				}



			}

		}

		if (strcmp(line,"#data")==0)
		{
			data_found=TRUE;
		}

	}while(!feof(file));
	fclose(file);
}
