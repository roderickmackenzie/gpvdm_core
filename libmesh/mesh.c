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
#include <mesh.h>
#include "inp.h"
#include "util.h"
#include "const.h"
#include "hard_limit.h"
#include <log.h>
#include <cal_path.h>
#include <lang.h>
#include <shape.h>

void mesh_cpy(struct simulation *sim,struct mesh *out,struct mesh *in)
{
	int i;
	int ii;

	mesh_free(out);

	out->nlayers= in->nlayers;
	out->remesh= in->remesh;


	out->layers = malloc (out->nlayers * sizeof(struct mesh_layer));

	for (i=0;i<out->nlayers;i++)
	{
		out->layers[i].len=in->layers[i].len;
		out->layers[i].n_points=in->layers[i].n_points;
		out->layers[i].mul=in->layers[i].mul;
		out->layers[i].left_right=in->layers[i].left_right;
	}


	for (i=0;i<out->nlayers;i++)
	{
		out->layers[i].n_points=in->layers[i].n_points;

		in->layers[i].dmesh=malloc ( out->layers[i].n_points * sizeof(long double));

		for (ii=0;ii>out->layers[i].n_points;ii++)
		{
			out->layers[i].dmesh[ii]=in->layers[i].dmesh[ii];
		}

	}

	out->tot_points=in->tot_points;

}

void mesh_check_y(struct simulation *sim,struct mesh *in,struct device *dev)
{
int y=0;
gdouble mesh_len=0.0;


	for (y=0;y<in->nlayers;y++)
	{
		mesh_len+=in->layers[y].len;
	}

	if (fabs(dev->ylen-mesh_len)>1e-14)
	{
		printf_log(sim,"calling remesh %Le %Le\n",dev->ylen,mesh_len);
		//getchar();
		mesh_remesh(sim,in,dev);
		printf_log(sim,"Warning: Length of epitaxy and computational mesh did not match, so I remesshed the device.\n");
	}
}

void mesh_remesh(struct simulation *sim,struct mesh *in,struct device *dev)
{
	char device_file_path[1000];

	if (in->remesh==TRUE)
	{

		in->nlayers=1;
		in->layers[0].len=dev->ylen;
		if (dev->my_epitaxy.electrical_layers==1)
		{
			in->layers[0].n_points=10;
		}else
		{
			in->layers[0].n_points=40;
		}
		in->layers[0].mul=1.0;
		in->layers[0].left_right=FALSE;
		mesh_save(sim,"mesh_y.inp",in);

		device_free(sim,dev);
		mesh_obj_free(sim,&(dev->mesh_data));

		join_path(2,device_file_path,get_input_path(sim),"epitaxy.inp");
		epitaxy_free(sim,&(dev->my_epitaxy));
		epitaxy_load(sim,&(dev->my_epitaxy),device_file_path);

		mesh_obj_load(sim,&(dev->mesh_data));

		dev->ns.dim.zmeshpoints=dev->mesh_data.meshdata_z.tot_points;
		dev->ns.dim.xmeshpoints=dev->mesh_data.meshdata_x.tot_points;
		dev->ns.dim.ymeshpoints=dev->mesh_data.meshdata_y.tot_points;

		device_get_memory(sim,dev);
		mesh_build(sim,dev);
		//printf("remesh\n");
		//getchar();
	}else
	{
		ewe(sim,"%s\n",_("The mesh does not match the device length and I am not alowed to remesh it"));
	}


}

void mesh_save(struct simulation *sim,char *file_name,struct mesh *in)
{
	int i=0;
	char buffer[2000];
	char temp[2000];
	char full_file_name[200];

	strcpy(buffer,"");
	strcat(buffer,"#remesh_enable\n");
	strcat(buffer,"True\n");
	strcat(buffer,"#mesh_layers\n");

	sprintf(temp,"%d\n",in->nlayers);
	strcat(buffer,temp);

	for (i=0;i<in->nlayers;i++)
	{
		strcat(buffer,"#mesh_layer_length0\n");

		sprintf(temp,"%Le\n",in->layers[i].len);
		strcat(buffer,temp);

		strcat(buffer,"#mesh_layer_points0\n");

		sprintf(temp,"%d\n",(int)(in->layers[i].n_points));
		strcat(buffer,temp);

		strcat(buffer,"#mesh_layer_mul0\n");

		sprintf(temp,"%d\n",(int)(in->layers[i].mul));
		strcat(buffer,temp);

		strcat(buffer,"#mesh_layer_left_right0\n");

		sprintf(temp,"%d\n",(int)(in->layers[i].left_right));
		strcat(buffer,temp);
	}

	strcat(buffer,"#ver\n");
	strcat(buffer,"1.0\n");
	strcat(buffer,"#end\n");

	join_path(2,full_file_name,get_input_path(sim),file_name);
	printf_log(sim,"Write new mesh to: %s\n",full_file_name);
	zip_write_buffer(sim,full_file_name,buffer, strlen(buffer));

}

void mesh_free(struct mesh *in)
{
	int i=0;

	if (in->nlayers==-1)
	{
		for (i=0;i<in->nlayers;i++)
		{
			free(in->layers[i].dmesh);
		}

		free(in->layers);
		mesh_init(in);
	}


}

void mesh_init(struct mesh *in)
{
	in->layers= NULL;
	in->nlayers=-1;
	in->remesh=-1;
}


void mesh_malloc_sub_mesh(struct simulation * sim, struct mesh *in)
{
	int i;
	long double pos;

	pos=0.0;
	int points=0;
	long double dx=0.0;

	for (i=0;i<in->nlayers;i++)
	{
		dx=in->layers[i].len/((long double)in->layers[i].n_points);
		pos=0.0;
		in->layers[i].n_points=0;
		while(pos<in->layers[i].len)
		{
			pos+=dx/2.0;
			//printf("%Le %Le\n",pos,(*meshdata)[i].len);
			//getchar();
			if (pos>in->layers[i].len)
			{
				break;
			}

			in->layers[i].n_points++;
			points++;
			pos+=dx/2.0;
			dx=dx*in->layers[i].mul;
		}

		in->layers[i].dmesh=malloc ( in->layers[i].n_points * sizeof(long double));

	}

	in->tot_points=points;

}

void mesh_gen_simple(struct simulation * sim, struct mesh *in,long double len,int points)
{
	int i;

	if (points<=0)
	{
		ewe(sim,"%s\n",_("Can't generate a mesh with zero points"));
	}
	in->remesh=FALSE;
	in->nlayers=1;

	in->layers = malloc (in->nlayers * sizeof(struct mesh_layer));

	i=0;
	in->layers[i].len=len;
	in->layers[i].n_points=points;
	in->layers[i].dx=in->layers[i].len/((long double)in->layers[i].n_points);
	in->layers[i].mul=1.0;
	in->layers[i].left_right=TRUE;


	mesh_malloc_sub_mesh(sim, in);

}

void mesh_gen_graded(struct simulation * sim, struct mesh *in,long double len,int points)
{
	int i;

	in->remesh=FALSE;
	in->nlayers=2;
	len/=2.0;

	in->layers = malloc (in->nlayers * sizeof(struct mesh_layer));

	i=0;
	in->layers[i].len=len;
	in->layers[i].n_points=points;
	in->layers[i].dx=in->layers[i].len/((long double)in->layers[i].n_points);
	in->layers[i].mul=1.1;
	in->layers[i].left_right=FALSE;

	i++;
	in->layers[i].len=len;
	in->layers[i].n_points=points;
	in->layers[i].dx=in->layers[i].len/((long double)in->layers[i].n_points);
	in->layers[i].mul=1.1;
	in->layers[i].left_right=TRUE;

	mesh_malloc_sub_mesh(sim, in);

}
void mesh_load_file(struct simulation * sim, struct mesh *in,char *file)
{
	int i;
	struct inp_file inp;
	char token0[200];
	char token1[200];
	//char token1[200];
	char val[200];
	//long double dx=0.0;
	//int points=0;
	//long double pos=0.0;

	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_input_path(sim),file);
	inp_check(sim,&inp,1.0);
	inp_reset_read(sim,&inp);

	inp_get_string(sim,&inp);			//remesh
	strcpy(val,inp_get_string(sim,&inp));
	in->remesh=english_to_bin(sim,val);

	inp_get_string(sim,&inp);			//layers
	sscanf(inp_get_string(sim,&inp),"%d",&(in->nlayers));

	in->layers = malloc (in->nlayers * sizeof(struct mesh_layer));

	for (i=0;i<in->nlayers;i++)
	{
		sscanf(inp_get_string(sim,&inp),"%s",token0);
		sscanf(inp_get_string(sim,&inp),"%Lf",&(in->layers[i].len));

		sscanf(inp_get_string(sim,&inp),"%s",token1);
		sscanf(inp_get_string(sim,&inp),"%d",&(in->layers[i].n_points));
		in->layers[i].dx=in->layers[i].len/((long double)in->layers[i].n_points);

		sscanf(inp_get_string(sim,&inp),"%s",token1);
		sscanf(inp_get_string(sim,&inp),"%Lf",&(in->layers[i].mul));

		sscanf(inp_get_string(sim,&inp),"%s",token1);
		sscanf(inp_get_string(sim,&inp),"%s",val);
		in->layers[i].left_right=english_to_bin(sim,val);

		in->layers[i].len=fabs(in->layers[i].len);
		hard_limit(sim,token0,&(in->layers[i].len));

	}

	inp_free(sim,&inp);

	mesh_malloc_sub_mesh(sim, in);

}

long double mesh_to_dim(struct simulation *sim,struct dimensions *dim, struct mesh *in,char xyz)
{

	long double *mesh;
	long double *dmesh;
	long double ret_len=0.0;

	dim_free_xyz(dim,xyz);

	if (xyz=='x')
	{
		dim->xmeshpoints=in->tot_points;
		dim_alloc_xyz(dim,'x');
		mesh=dim->xmesh;
		dmesh=dim->dxmesh;
	}else
	if (xyz=='y')
	{
		dim->ymeshpoints=in->tot_points;
		dim_alloc_xyz(dim,'y');
		mesh=dim->ymesh;
		dmesh=dim->dymesh;
	}else
	if (xyz=='z')
	{
		dim->zmeshpoints=in->tot_points;
		dim_alloc_xyz(dim,'z');
		mesh=dim->zmesh;
		dmesh=dim->dzmesh;
	}

	int pos=0;
	int i=0;
	int ii=0;
	//int z=0;
	//int x=0;
	gdouble dpos=0.0;
	long double dx=0.0;
	long double len=0.0;

	len=0.0;
	for (i=0;i<in->nlayers;i++)
	{
		pos=0;
		dpos=0.0;
		dx=in->layers[i].dx;
		//printf("going to build %ld\n",meshdata[i].n_points);
		//getchar();
		for (ii=0;ii<in->layers[i].n_points;ii++)
		{
			dpos+=dx/2.0;
			in->layers[i].dmesh[ii]=dpos;
			dpos+=dx/2.0;
			dx*=in->layers[i].mul;
			pos++;
		}
		len+=in->layers[i].len;
	}

	ret_len=len;

	len=0.0;
	pos=0;
	for (i=0;i<in->nlayers;i++)
	{
		for (ii=0;ii<in->layers[i].n_points;ii++)
		{
			//printf("%d %d %d %Le\n",pos,in->layers[i].n_points,in->nlayers,in->layers[i].len);
			if (in->layers[i].left_right==FALSE)
			{
				mesh[pos]=len+in->layers[i].dmesh[ii];
			}else
			{
				mesh[pos]=len+in->layers[i].len-in->layers[i].dmesh[in->layers[i].n_points-1-ii];
			}
			//printf("%c %ld %Le %d %d %ld\n",direction,pos,mesh[pos],i,ii,meshdata[i].n_points);
			pos++;
		}
		len+=in->layers[i].len;
	}

	//build dmesh
	long double last=0.0;
	long double next=0.0;

	for (i=0;i<in->tot_points;i++)
	{
		if ((in->tot_points-1)==i)
		{
			next=ret_len;
		}else
		{
			next=(mesh[i]+mesh[i+1])/2.0;
		}

		dx=next-last;
		dmesh[i]=dx;
		last=next;
	}


	return ret_len;

}

void mesh_dump(struct simulation *sim,struct dimensions *dim)
{
	int x=0;

	for (x=0;x<dim->xmeshpoints;x++)
	{
		printf("%Le\n",dim->xmesh[x]);
	}
}

void mesh_dump_y(struct simulation *sim,struct dimensions *dim)
{
	int y=0;

	for (y=0;y<dim->ymeshpoints;y++)
	{
		printf("%Le\n",dim->ymesh[y]);
	}
}


void mesh_build(struct simulation *sim,struct device *in)
{

	in->zlen=mesh_to_dim(sim, &(in->ns.dim), &(in->mesh_data.meshdata_z),'z');
	in->xlen=mesh_to_dim(sim, &(in->ns.dim), &(in->mesh_data.meshdata_x),'x');
	in->ylen=mesh_to_dim(sim, &(in->ns.dim), &(in->mesh_data.meshdata_y),'y');
	//dim_cpy(&(in->dim_max),&(in->ns.dim));
}

void mesh_numerate_points(struct simulation *sim,struct device *in)
{

	int shape=0;
	int z=0;
	int x=0;
	int y=0;

	struct dimensions *dim=&(in->ns.dim);

	gdouble dpos=0.0;

	//len=0.0;
	for (y=0;y<dim->ymeshpoints;y++)
	{
		dpos=dim->ymesh[y];
		in->imat[0][0][y]=epitaxy_get_electrical_material_layer(&(in->my_epitaxy),dpos);
		in->imat_epitaxy[0][0][y]=epitaxy_get_epitaxy_layer_using_electrical_pos(&(in->my_epitaxy),dpos);

		for (z=0;z<dim->zmeshpoints;z++)
		{
			for (x=0;x<dim->xmeshpoints;x++)
			{
				in->imat[z][x][y]=in->imat[0][0][y];
				in->imat_epitaxy[z][x][y]=in->imat_epitaxy[0][0][y];
				shape=shape_get_index(sim,&(in->my_epitaxy),dim->xmesh[x],dim->ymesh[y],dim->zmesh[z]);
				if (shape!=-1)
				{
					in->imat[z][x][y]=shape;
				}
				//printf("%d\n",shape);
			}
		}

	}

}
void mesh_cal_layer_widths(struct device *in)
{
int i;
int cur_i=in->imat[0][0][0];

in->layer_start[cur_i]=0.0;
struct dimensions *dim=&in->ns.dim;

for (i=0;i<dim->ymeshpoints;i++)
{
	if ((in->imat[0][0][i]!=cur_i)||(i==(dim->ymeshpoints-1)))
	{
		in->layer_stop[cur_i]=dim->ymesh[i-1];//+(dim->ymesh[i]-dim->ymesh[i-1])/2;
		if (i==(dim->ymeshpoints-1))
		{
			break;
		}
		cur_i=in->imat[0][0][i];
		in->layer_start[cur_i]=dim->ymesh[i];//-(dim->ymesh[i]-dim->ymesh[i-1])/2;
	}
//printf_log("%d\n",in->imat[i]);
}
}
