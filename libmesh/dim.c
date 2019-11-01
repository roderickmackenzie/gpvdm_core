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

/** @file dim.c
@brief Dimension object
*/

#include <device.h>
#include "mesh.h"
#include "inp.h"
#include "util.h"
#include "const.h"
#include "hard_limit.h"
#include <log.h>
#include <cal_path.h>
#include <lang.h>
#include <shape.h>


void dim_init_xyz(struct dimensions *dim,char xyz)
{
	if (xyz=='x')
	{
		dim->xmesh= NULL;
		dim->dxmesh= NULL;
		dim->xmeshpoints=-1;
	}else
	if (xyz=='y')
	{
		dim->ymesh= NULL;
		dim->dymesh= NULL;
		dim->ymeshpoints=-1;
	}else
	if (xyz=='z')
	{
		dim->zmesh= NULL;
		dim->dzmesh= NULL;
		dim->zmeshpoints=-1;
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
			free(dim->dxmesh);
			dim_init_xyz(dim,'x');
		}
	}else
	if (xyz=='y')
	{
		if (dim->ymesh!=NULL)
		{
			free(dim->ymesh);
			free(dim->dymesh);
			dim_init_xyz(dim,'y');
		}
	}else
	if (xyz=='z')
	{
		if (dim->zmesh!=NULL)
		{
			free(dim->zmesh);
			free(dim->dzmesh);
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
		dim->xmesh = (gdouble *) malloc(dim->xmeshpoints * sizeof(gdouble));
		memset(dim->xmesh, 0, dim->xmeshpoints * sizeof(gdouble));

		dim->dxmesh = (gdouble *) malloc(dim->xmeshpoints * sizeof(gdouble));
		memset(dim->dxmesh, 0, dim->xmeshpoints * sizeof(gdouble));
	}else
	if (xyz=='y')
	{
		dim->ymesh = (gdouble *) malloc(dim->ymeshpoints * sizeof(gdouble));
		memset(dim->ymesh, 0, dim->ymeshpoints * sizeof(gdouble));

		dim->dymesh = (gdouble *) malloc(dim->ymeshpoints * sizeof(gdouble));
		memset(dim->dymesh, 0, dim->ymeshpoints * sizeof(gdouble));
	}else
	if (xyz=='z')
	{
		dim->zmesh = (gdouble *) malloc(dim->zmeshpoints * sizeof(gdouble));
		memset(dim->zmesh, 0, dim->zmeshpoints * sizeof(gdouble));

		dim->dzmesh = (gdouble *) malloc(dim->zmeshpoints * sizeof(gdouble));
		memset(dim->dzmesh, 0, dim->zmeshpoints * sizeof(gdouble));
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
	out->xmeshpoints=in->xmeshpoints;
	out->ymeshpoints=in->ymeshpoints;
	out->zmeshpoints=in->zmeshpoints;
	dim_alloc(out);

	for (x=0;x<out->xmeshpoints;x++)
	{
		out->xmesh[x]=in->xmesh[x];
		out->dxmesh[x]=in->dxmesh[x];
	}

	for (y=0;y<out->ymeshpoints;y++)
	{
		out->ymesh[y]=in->ymesh[y];
		out->dymesh[y]=in->dymesh[y];
	}

	for (z=0;z<out->zmeshpoints;z++)
	{
		out->zmesh[z]=in->zmesh[z];
		out->dzmesh[z]=in->dzmesh[z];
	}

}

