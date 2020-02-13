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

/** @file dump.c
@brief go and dump stuff, what is dumped depends on which options have been set
*/

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <exp.h>
#include "sim.h"
#include "dump.h"
#include <cal_path.h>
#include <pl.h>
#include <probe.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include <inp.h>
#include <color.h>
#include <memory.h>
#include <ray_fun.h>
#include <device_fun.h>


static int unused __attribute__((unused));


void dump_clean_cache_files(struct simulation* sim)
{
struct inp_file inp;
char temp[200];
char cach_dir[PATH_MAX];


	if (get_dump_status(sim,dump_remove_dos_cache)==TRUE)
	{
		join_path(2, cach_dir,get_input_path(sim),"cache");
		remove_dir(sim,cach_dir);
	}

	inp_init(sim,&inp);
	if (inp_load(sim, &inp , "delete_files.inp")==0)
	{

		inp_reset_read(sim,&inp);
		strcpy(temp,inp_get_string(sim,&inp));
		if (strcmp(temp,"#begin")!=0)
		{
			return;
		}

		while(1)
		{
			strcpy(temp,inp_get_string(sim,&inp));
			if (strcmp(temp,"#end")==0)
			{
				break;
			}else
			{
				remove(temp);
			}

		}

		inp_free(sim,&inp);
	}



}

void dump_init(struct simulation *sim,struct device* in)
{
in->snapshot_number=0;
set_dump_status(sim,dump_lock, FALSE);
}

void buffer_add_3d_to_2d_projection(struct simulation *sim,struct dat_file *buf,struct device *in,gdouble ***data)
{
int x=0;
int y=0;
int z=0;

gdouble xpos=0.0;
gdouble ypos=0.0;
gdouble zpos=0.0;
long double tot=0.0;

struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

char string[200];
if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#data\n");
	buffer_add_string(buf,string);
}

if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			tot=0.0;
			for (y=0;y<dim->ylen;y++)
			{
				tot+=data[z][x][y];
			}

			sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],dim->zmesh[z],tot/((long double)dim->ylen));
			buffer_add_string(buf,string);
		}
	}
}else
if ((dim->xlen>1)&&(dim->ylen>1))
{
	z=0;
	for (x=0;x<dim->xlen;x++)
	{
		tot=0.0;
		for (y=0;y<dim->ylen;y++)
		{
			tot+=data[z][x][y];
		}
		sprintf(string,"%Le %Le\n",dim->xmesh[x],tot/((long double)dim->ylen));
		buffer_add_string(buf,string);
	}
}else
{
	x=0;
	z=0;
	tot=0.0;
	for (y=0;y<dim->ylen;y++)
	{
		tot+=data[z][x][y];
	}

	sprintf(string,"%Le %Le\n",dim->xmesh[x],tot/((long double)dim->ylen));
	buffer_add_string(buf,string);
}

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#end\n");
	buffer_add_string(buf,string);
}

}

void buffer_add_3d_data(struct simulation *sim,struct dat_file *buf,struct dimensions *dim,gdouble ***data)
{
int x=0;
int y=0;
int z=0;

gdouble xpos=0.0;
gdouble ypos=0.0;
gdouble zpos=0.0;

char string[200];
if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#data\n");
	buffer_add_string(buf,string);
}

if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				sprintf(string,"%Le %Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],dim->zmesh[z],data[z][x][y]);
				buffer_add_string(buf,string);
			}
		}
	}
}else
if ((dim->xlen>1)&&(dim->ylen>1))
{
	z=0;
	for (x=0;x<dim->xlen;x++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],data[z][x][y]);
			buffer_add_string(buf,string);
		}
		buffer_add_string(buf,"\n");
	}
}else
{
	x=0;
	z=0;
	for (y=0;y<dim->ylen;y++)
	{
		sprintf(string,"%Le %Le\n",dim->ymesh[y],data[z][x][y]);
		buffer_add_string(buf,string);
	}
}

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#end\n");
	buffer_add_string(buf,string);
}

}


void buffer_add_zxy_rgb_data(struct simulation *sim,struct dat_file *buf,struct dimensions *dim,gdouble ***data)
{
int x=0;
int z=0;
int y=0;

gdouble xpos=0.0;
gdouble zpos=0.0;
long double y_tot=0;

long double X;
long double Y;
long double Z;
int R;
int G;
int B;

struct istruct luminescence_tot;

char string[200];
struct dimensions XYZ_dim;
long double ***XYZ;
long double max;
dim_init(&XYZ_dim);
dim_cpy(&XYZ_dim,dim);
dim_free_xyz(&XYZ_dim,'y');
XYZ_dim.ylen=3;
dim_alloc_xyz(&XYZ_dim,'y');
XYZ_dim.ymesh[0]=0.0;
XYZ_dim.ymesh[1]=1.0;
XYZ_dim.ymesh[2]=2.0;

malloc_zxy_gdouble(&(XYZ_dim),&XYZ);


	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#data\n");
		buffer_add_string(buf,string);
	}

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{

			inter_init(sim,&luminescence_tot);

			for (y=0;y<dim->ylen;y++)
			{
				inter_append(&luminescence_tot,dim->ymesh[y],data[z][x][y]);
			}

			color_cie_cal_XYZ(sim,&X,&Y,&Z,&luminescence_tot,FALSE);
			XYZ[z][x][0]=X;
			XYZ[z][x][1]=Y;
			XYZ[z][x][2]=Z;

			inter_free(&luminescence_tot);
		}
	}

	max=zx_y_max_gdouble(&XYZ_dim, XYZ,1);
	zxy_div_gdouble(&XYZ_dim, XYZ, max);

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{

			X=XYZ[z][x][0];
			Y=XYZ[z][x][1];
			Z=XYZ[z][x][2];

			color_XYZ_to_rgb(&R,&G,&B,X,Y,Z);

			R*=Y;
			G*=Y;
			B*=Y;

			sprintf(string,"%Le %Le %.2x%.2x%.2x\n",dim->zmesh[z],dim->xmesh[x],R,G,B);
			buffer_add_string(buf,string);

		}

		sprintf(string,"\n");
		buffer_add_string(buf,string);
	}

	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#end\n");
		buffer_add_string(buf,string);
	}

zx_y_quick_dump("XYZ", XYZ, &(XYZ_dim));

free_zxy_gdouble(&XYZ_dim,&XYZ);
dim_free(&XYZ_dim);
}

void buffer_add_yl_light_data(struct simulation *sim,struct dat_file *buf,struct dim_light *dim,long double ****data,long double shift, int z, int x)
{
int y=0;
int l=0;

long double xpos=0.0;
long double zpos=0.0;

char string[200];

	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#data\n");
		buffer_add_string(buf,string);
	}

	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(string,"%Le %Le %Le\n",dim->l[l],dim->y[y]-shift,data[z][x][y][l]);

			buffer_add_string(buf,string);
		}

		buffer_add_string(buf,"\n");
	}

	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#end\n");
		buffer_add_string(buf,string);
	}

}

void buffer_add_zx_data(struct simulation *sim,struct dat_file *buf,struct dimensions *dim,gdouble **data)
{
int x=0;
int z=0;

gdouble xpos=0.0;
gdouble zpos=0.0;

char string[200];

	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#data\n");
		buffer_add_string(buf,string);
	}


	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			sprintf(string,"%Le %Le %Le \n",dim->zmesh[z],dim->xmesh[x],data[z][x]);
			buffer_add_string(buf,string);
		}
	}

	if (get_dump_status(sim,dump_write_headers)==TRUE)
	{
		sprintf(string,"#end\n");
		buffer_add_string(buf,string);
	}

}

void buffer_add_3d_device_data_int(struct simulation *sim,struct dat_file *buf,struct device *in,int ***data)
{
int x=0;
int y=0;
int z=0;

gdouble xpos=0.0;
gdouble ypos=0.0;
gdouble zpos=0.0;

struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

char string[200];
if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#data\n");
	buffer_add_string(buf,string);
}

if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				sprintf(string,"%Le %Le %Le %d\n",dim->xmesh[x],dim->ymesh[y],dim->zmesh[z],data[z][x][y]);
				buffer_add_string(buf,string);
			}
		}
	}
}else
if ((dim->xlen>1)&&(dim->ylen>1))
{
	z=0;
	for (x=0;x<dim->xlen;x++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(string,"%Le %Le %d\n",dim->xmesh[x],dim->ymesh[y],data[z][x][y]);
			buffer_add_string(buf,string);
		}
		buffer_add_string(buf,"\n");
	}
}else
{
	x=0;
	z=0;
	for (y=0;y<dim->ylen;y++)
	{
		sprintf(string,"%Le %d\n",dim->ymesh[y],data[z][x][y]);
		buffer_add_string(buf,string);
	}
}

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#end\n");
	buffer_add_string(buf,string);
}

}

void buffer_add_2d_device_data_int(struct simulation *sim,struct dat_file *buf,struct device *in,int **data)
{
int x=0;
int z=0;
struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

gdouble xpos=0.0;
gdouble zpos=0.0;

char string[200];
if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#data\n");
	buffer_add_string(buf,string);
}

if ((dim->xlen>1)&&(dim->zlen>1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
				sprintf(string,"%Le %Le %d\n",dim->xmesh[x],dim->zmesh[z],data[z][x]);
				buffer_add_string(buf,string);
		}

		buffer_add_string(buf,"\n");
	}
}else
if ((dim->xlen>1))
{
	z=0;
	for (x=0;x<dim->xlen;x++)
	{
		sprintf(string,"%Le %d\n",dim->xmesh[x],data[z][x]);
		buffer_add_string(buf,string);
	}
}else
{
	sprintf(string,"%Le %d\n",dim->ymesh[0],data[0][0]);
	buffer_add_string(buf,string);
}

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#end\n");
	buffer_add_string(buf,string);
}

}

void buffer_add_3d_device_data_including_boundaries(struct simulation *sim,struct dat_file *buf,struct device *in,gdouble ***data,long double **left,long double **right)
{
int x=0;
int y=0;
int z=0;

gdouble xpos=0.0;
gdouble ypos=0.0;
gdouble zpos=0.0;

char string[200];
struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#data\n");
	buffer_add_string(buf,string);
}

if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],(long double)0.0,left[z][x]);
			buffer_add_string(buf,string);

			for (y=0;y<dim->ylen;y++)
			{
				sprintf(string,"%Le %Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],dim->zmesh[z],data[z][x][y]);
				buffer_add_string(buf,string);
			}

			sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],in->ylen,right[z][x]);
			buffer_add_string(buf,string);

		}
	}
}else
if ((dim->xlen>1)&&(dim->ylen>1))
{
	z=0;
	for (x=0;x<dim->xlen;x++)
	{
		sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],(long double)0.0,left[z][x]);
		buffer_add_string(buf,string);

		for (y=0;y<dim->ylen;y++)
		{
			sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],dim->ymesh[y],data[z][x][y]);
			buffer_add_string(buf,string);
		}

		sprintf(string,"%Le %Le %Le\n",dim->xmesh[x],in->ylen,right[z][x]);
		buffer_add_string(buf,string);

		buffer_add_string(buf,"\n");
	}
}else
{
	x=0;
	z=0;
	sprintf(string,"%Le %Le\n",(long double)0.0,left[z][x]);
	buffer_add_string(buf,string);

	for (y=0;y<dim->ylen;y++)
	{
		sprintf(string,"%Le %Le\n",dim->ymesh[y],data[z][x][y]);
		buffer_add_string(buf,string);
	}

	sprintf(string,"%Le %Le\n",in->ylen,right[z][x]);
	buffer_add_string(buf,string);

}

if (get_dump_status(sim,dump_write_headers)==TRUE)
{
	sprintf(string,"#end\n");
	buffer_add_string(buf,string);
}

}

void buffer_set_graph_type(struct dat_file *buf,struct device *in)
{
	struct dimensions *dim=&in->ns.dim;

	if ((dim->xlen>1)&&(dim->ylen>1)&&(dim->zlen>1))
	{
		strcpy(buf->type,"4d");
	}else
	if ((dim->xlen>1)&&(dim->ylen>1))
	{
		strcpy(buf->type,"3d");
	}else
	{
		strcpy(buf->type,"xy");
	}
}

void dump_write_to_disk(struct simulation *sim,struct device* in)
{
char temp[200];
char postfix[100];
char out_dir[PATH_MAX];
char sim_name[PATH_MAX];
struct dat_file buf;
buffer_init(&buf);

strextract_name(sim_name,in->simmode);


int dumped=FALSE;
FILE* out;

	sprintf(postfix,"%d",in->snapshot_number);
	//if ((get_dump_status(sim,dump_pl)==TRUE)||(get_dump_status(sim,dump_energy_slice_switch)==TRUE)||(get_dump_status(sim,dump_1d_slices)==TRUE)||(get_dump_status(sim,dump_optical_probe_spectrum)==TRUE))
	//{

	dump_make_snapshot_dir(sim,out_dir ,in->time, get_equiv_V(sim,in), in->fx, in->snapshot_number);
	//}


	if ((in->dump_1d_slice_zpos!=-1)&&(in->dump_1d_slice_xpos!=-1))
	{
		dump_device_map(sim,out_dir,in);
		dumped=TRUE;
	}


	if (get_dump_status(sim,dump_1d_slices)==TRUE)
	{
		dump_1d_slice(sim,in,out_dir);
		dumped=TRUE;
	}

	if (get_dump_status(sim,dump_energy_slice_switch)==TRUE)
	{
		dump_energy_slice(sim,out_dir,in);
		dumped=TRUE;
	}



	if (dumped==TRUE)
	{
		in->snapshot_number++;
	}


}


