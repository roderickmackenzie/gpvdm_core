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


void dim_init(struct dimensions *dim)
{
	dim->ymesh= NULL;
	dim->xmesh= NULL;
	dim->zmesh= NULL;

	dim->dymesh=NULL;
	dim->dxmesh=NULL;
	dim->dzmesh=NULL;

	dim->zmeshpoints=-1;
	dim->xmeshpoints=-1;
	dim->ymeshpoints=-1;
	dim->srh_bands=-1;
}

void dim_free(struct dimensions *dim)
{
	if (dim->ymesh!=NULL)
	{
		free(dim->ymesh);
		free(dim->xmesh);
		free(dim->zmesh);

		free(dim->dymesh);
		free(dim->dxmesh);
		free(dim->dzmesh);

		dim_init(dim);
	}
}

void dim_alloc(struct dimensions *dim)
{

	dim->zmesh = (gdouble *) malloc(dim->zmeshpoints * sizeof(gdouble));
	memset(dim->zmesh, 0, dim->zmeshpoints * sizeof(gdouble));

	dim->xmesh = (gdouble *) malloc(dim->xmeshpoints * sizeof(gdouble));
	memset(dim->xmesh, 0, dim->xmeshpoints * sizeof(gdouble));

	dim->ymesh = (gdouble *) malloc(dim->ymeshpoints * sizeof(gdouble));
	memset(dim->ymesh, 0, dim->ymeshpoints * sizeof(gdouble));

	dim->dzmesh = (gdouble *) malloc(dim->zmeshpoints * sizeof(gdouble));
	memset(dim->dzmesh, 0, dim->zmeshpoints * sizeof(gdouble));

	dim->dxmesh = (gdouble *) malloc(dim->xmeshpoints * sizeof(gdouble));
	memset(dim->dxmesh, 0, dim->xmeshpoints * sizeof(gdouble));

	dim->dymesh = (gdouble *) malloc(dim->ymeshpoints * sizeof(gdouble));
	memset(dim->dymesh, 0, dim->ymeshpoints * sizeof(gdouble));

	printf("%d %d %d %d %d %d\n",dim->ymesh,dim->xmesh,dim->zmesh,dim->dymesh,dim->dxmesh,dim->dzmesh);

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

