//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/** @file dim.c
@brief Dimension object
*/

#include <device.h>
#include "mesh.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include "hard_limit.h"
#include <log.h>
#include <cal_path.h>
#include <lang.h>
#include <shape.h>
#include <device_fun.h>
#include <math.h>

void dim_init_xyz(struct dimensions *dim,char xyz)
{
	if (xyz=='x')
	{
		dim->xmesh= NULL;
		dim->dx= NULL;
		dim->xlen=-1;
	}else
	if (xyz=='y')
	{
		dim->ymesh= NULL;
		dim->dy= NULL;
		dim->ylen=-1;
	}else
	if (xyz=='z')
	{
		dim->zmesh= NULL;
		dim->dz= NULL;
		dim->zlen=-1;
	}

}

void dim_init(struct dimensions *dim)
{
	dim_init_xyz(dim,'x');
	dim_init_xyz(dim,'y');
	dim_init_xyz(dim,'z');
	dim->srh_bands=-1;
}

void dim_free_xyz(struct dimensions *dim,char xyz)
{
	if (xyz=='x')
	{
		if (dim->xmesh!=NULL)
		{
			free(dim->xmesh);
			free(dim->dx);
			dim_init_xyz(dim,'x');
		}
	}else
	if (xyz=='y')
	{
		if (dim->ymesh!=NULL)
		{
			free(dim->ymesh);
			free(dim->dy);
			dim_init_xyz(dim,'y');
		}
	}else
	if (xyz=='z')
	{
		if (dim->zmesh!=NULL)
		{
			free(dim->zmesh);
			free(dim->dz);
			dim_init_xyz(dim,'z');
		}
	}

}

void dim_free(struct dimensions *dim)
{
	dim_free_xyz(dim,'x');
	dim_free_xyz(dim,'y');
	dim_free_xyz(dim,'z');
	dim_init(dim);
}

void dim_alloc_xyz(struct dimensions *dim,char xyz)
{

	if (xyz=='x')
	{
		if (dim->xlen!=0)
		{
			dim->xmesh = (gdouble *) malloc(dim->xlen * sizeof(gdouble));
			memset(dim->xmesh, 0, dim->xlen * sizeof(gdouble));

			dim->dx = (gdouble *) malloc(dim->xlen * sizeof(gdouble));
			memset(dim->dx, 0, dim->xlen * sizeof(gdouble));
		}
	}else
	if (xyz=='y')
	{
		if (dim->ylen!=0)
		{
			dim->ymesh = (gdouble *) malloc(dim->ylen * sizeof(gdouble));
			memset(dim->ymesh, 0, dim->ylen * sizeof(gdouble));

			dim->dy = (gdouble *) malloc(dim->ylen * sizeof(gdouble));
			memset(dim->dy, 0, dim->ylen * sizeof(gdouble));
		}
	}else
	if (xyz=='z')
	{
		if (dim->zlen!=0)
		{
			dim->zmesh = (gdouble *) malloc(dim->zlen * sizeof(gdouble));
			memset(dim->zmesh, 0, dim->zlen * sizeof(gdouble));

			dim->dz = (gdouble *) malloc(dim->zlen * sizeof(gdouble));
			memset(dim->dz, 0, dim->zlen * sizeof(gdouble));
		}
	}


}

void dim_alloc_gen_untiy_mesh_x(struct dimensions *dim)
{

	int x=0;
	long double xpos=0.0;
	long double dx=1.0/dim->xlen;

	for (x=0;x<dim->xlen;x++)
	{
		dim->xmesh[x]=xpos;
		xpos+=dx;
	}
}

void dim_alloc_gen_untiy_mesh_z(struct dimensions *dim)
{
	int z=0;
	long double zpos=0.0;
	long double dz=1.0/dim->zlen;

	for (z=0;z<dim->zlen;z++)
	{
		dim->zmesh[z]=zpos;
		zpos+=dz;
	}
}

void dim_alloc(struct dimensions *dim)
{
	dim_alloc_xyz(dim,'x');
	dim_alloc_xyz(dim,'y');
	dim_alloc_xyz(dim,'z');

}

void dim_swap(struct dimensions *out,struct dimensions *in)
{
	struct dimensions tmp;
	dim_init(&tmp);
	dim_cpy(&tmp,out);
	dim_cpy(out,in);
	dim_cpy(in,&tmp);
	dim_free(&tmp);

}

void dim_cpy(struct dimensions *out,struct dimensions *in)
{

	int x=0;
	int y=0;
	int z=0;

	dim_free(out);
	out->xlen=in->xlen;
	out->ylen=in->ylen;
	out->zlen=in->zlen;

	dim_alloc(out);

	for (x=0;x<out->xlen;x++)
	{
		out->xmesh[x]=in->xmesh[x];
		out->dx[x]=in->dx[x];
	}

	for (y=0;y<out->ylen;y++)
	{
		out->ymesh[y]=in->ymesh[y];
		out->dy[y]=in->dy[y];
	}

	for (z=0;z<out->zlen;z++)
	{
		out->zmesh[z]=in->zmesh[z];
		out->dz[z]=in->dz[z];
	}

	out->srh_bands=in->srh_bands;
}

void dim_info_to_buf(struct dat_file *buf,struct dimensions *dim)
{
	long double mul_x=0.0;
	long double mul_y=0.0;
	long double mul_z=0.0;

	get_meter_dim(buf->x_units,&mul_x,dim->xmesh[dim->xlen-1]);
	get_meter_dim(buf->y_units,&mul_y,dim->ymesh[dim->ylen-1]);
	get_meter_dim(buf->z_units,&mul_z,dim->zmesh[dim->zlen-1]);
	buf->y_mul=mul_y;
	buf->x_mul=mul_x;
	buf->z_mul=mul_z;

	strcpy(buf->x_label,_("x-position"));
	strcpy(buf->y_label,_("y-position"));
	strcpy(buf->z_label,_("z-position"));

	buf->x=dim->xlen;
	buf->y=dim->ylen;
	buf->z=dim->zlen;

	buf->logscale_x=0;
	buf->logscale_y=0;

	if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
	{
		strcpy(buf->type,"zxy-d");
	}else
	if ((dim->xlen>1)&&(dim->ylen>1))
	{
		strcpy(buf->type,"3d");
	}else
	{
		strcpy(buf->type,"xy");
	}
}



long double dim_dl_xyz(struct dimensions *dim, int x0, int y0, int z0,int x1,int y1, int z1)
{
	long double pos_x0=dim->xmesh[x0];
	long double pos_y0=dim->ymesh[y0];
	long double pos_z0=dim->zmesh[z0];

	long double pos_x1=dim->xmesh[x1];
	long double pos_y1=dim->ymesh[y1];
	long double pos_z1=dim->zmesh[z1];

	return sqrtl(powl(pos_x0-pos_x1,2.0)+powl(pos_y0-pos_y1,2.0)+powl(pos_z0-pos_z1,2.0));
}


