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

/** @file mesh.c
@brief This builds the electrical mesh
*/

#include "device.h"
#include "mesh.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include "hard_limit.h"
#include <log.h>
#include <cal_path.h>
#include <lang.h>
#include <shape.h>


void mesh_obj_cpy(struct simulation *sim,struct mesh_obj *out,struct mesh_obj *in)
{
	mesh_cpy(sim,&(out->meshdata_x),&(in->meshdata_x));
	mesh_cpy(sim,&(out->meshdata_y),&(in->meshdata_y));
	mesh_cpy(sim,&(out->meshdata_z),&(in->meshdata_z));
}

void mesh_obj_free(struct simulation *sim,struct mesh_obj *in)
{
	mesh_free(&(in->meshdata_x));
	mesh_free(&(in->meshdata_y));
	mesh_free(&(in->meshdata_z));
}

void mesh_obj_init(struct mesh_obj *in)
{
	mesh_init(&(in->meshdata_x));
	mesh_init(&(in->meshdata_y));
	mesh_init(&(in->meshdata_z));
}


void mesh_obj_load(struct simulation *sim,struct mesh_obj *mesh)
{
	mesh_load_file(sim,&(mesh->meshdata_z),"mesh_z.inp");
	mesh_load_file(sim,&(mesh->meshdata_x),"mesh_x.inp");
	mesh_load_file(sim,&(mesh->meshdata_y),"mesh_y.inp");
}

