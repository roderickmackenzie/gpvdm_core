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

#include <stdio.h>
#include <ray.h>
#include <gpvdm_const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <dat_file.h>
#include <string.h>
#include <util.h>


/** @file triangles.c
	@brief Basic low level triangle functions
*/

void triangle_load_from_file(struct simulation *sim,struct triangles *in,char *file_name)
{
char line[10000];
struct dat_file dat;
dat_file_load_info(sim,&dat,file_name);
in->len=dat.y;
in->data=(struct triangle*)malloc(sizeof(struct triangle)*in->len);

FILE *file;
file=fopen(file_name,"r");
int point=0;
int items=0;
if (file == NULL)
{
	ewe(sim,"Can not load the file %s\n",file_name);
}
int data_start=FALSE;
do
{
	memset(line,0,10000);
	fgets(line, 10000, file);

	if (data_start==FALSE)
	{
		if (strcmp_begin(line,"#begin")==0)
		{
			data_start=TRUE;
		}
	}else
	{
		if (line[0]!='#')
		{
			if (point==0)
			{
				sscanf(line,"%le %le %le",&(in->data[items].xy0.z),&(in->data[items].xy0.x),&(in->data[items].xy0.y));

			}else
			if (point==1)
			{
				sscanf(line,"%le %le %le",&(in->data[items].xy1.z),&(in->data[items].xy1.x),&(in->data[items].xy1.y));
			}else
			if (point==2)
			{
				sscanf(line,"%le %le %le",&(in->data[items].xy2.z),&(in->data[items].xy2.x),&(in->data[items].xy2.y));
			}

			point++;

			if (point==5)
			{
				point=0;
				items++;
			}
		}


		if (items>in->len)
		{
			break;
		}

	}


}while(!feof(file));
fclose(file);

}


void triangles_free(struct triangles *tri)
{
	if (tri->data!=NULL)
	{
		free(tri->data);
	}

	tri->len=0;
	tri->max_len=0;

	tri->data=NULL;

}

void triangles_cpy(struct triangles *out,struct triangles *in)
{
	out->data=(struct triangle*)malloc(sizeof(struct triangle)*in->len);
	memcpy(out->data, in->data, in->len*sizeof(struct triangle));
	out->len=in->len;
}

void triangles_set_object_type(struct triangles *in,int object_type)
{
int i;
	for (i=0;i<in->len;i++)
	{
		in->data[i].object_type=object_type;
	}

}

void triangles_save(char *file_name,struct triangles *in)
{	
	int i;
	FILE *out;
	struct triangle *tri;
	out=fopen(file_name,"w");
	for (i=0;i<in->len;i++)
	{
		tri=&(in->data[i]);
		fprintf(out,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
		fprintf(out,"%le %le %le\n",tri->xy1.z,tri->xy1.x,tri->xy1.y);
		fprintf(out,"%le %le %le\n",tri->xy2.z,tri->xy2.x,tri->xy2.y);
		fprintf(out,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
		fprintf(out,"\n");
		fprintf(out,"\n");
	}

	fclose(out);
}

void triangles_find_min(struct vec *out,struct triangles *in)
{
	int i;
	double x=in->data[0].xy0.x;
	double y=in->data[0].xy0.y;
	double z=in->data[0].xy0.z;

	for (i=0;i<in->len;i++)
	{
		if (in->data[i].xy0.x<x)
		{
			x=in->data[i].xy0.x;
		}

		if (in->data[i].xy0.y<y)
		{
			y=in->data[i].xy0.y;
		}

		if (in->data[i].xy0.z<z)
		{
			z=in->data[i].xy0.z;
		}


		if (in->data[i].xy1.x<x)
		{
			x=in->data[i].xy1.x;
		}

		if (in->data[i].xy1.y<y)
		{
			y=in->data[i].xy1.y;
		}

		if (in->data[i].xy1.z<z)
		{
			z=in->data[i].xy1.z;
		}

		if (in->data[i].xy2.x<x)
		{
			x=in->data[i].xy2.x;
		}

		if (in->data[i].xy2.y<y)
		{
			y=in->data[i].xy2.y;
		}

		if (in->data[i].xy2.z<z)
		{
			z=in->data[i].xy2.z;
		}
	}

	out->x=x;
	out->y=y;
	out->z=z;

}

void triangles_find_max(struct vec *out,struct triangles *in)
{

	int i;

	double x=in->data[0].xy0.x;
	double y=in->data[0].xy0.y;
	double z=in->data[0].xy0.z;

	for (i=0;i<in->len;i++)
	{
		if (in->data[i].xy0.x>x)
		{
			x=in->data[i].xy0.x;
		}

		if (in->data[i].xy0.y>y)
		{
			y=in->data[i].xy0.y;
		}

		if (in->data[i].xy0.z>z)
		{
			z=in->data[i].xy0.z;
		}


		if (in->data[i].xy1.x>x)
		{
			x=in->data[i].xy1.x;
		}

		if (in->data[i].xy1.y>y)
		{
			y=in->data[i].xy1.y;
		}

		if (in->data[i].xy1.z>z)
		{
			z=in->data[i].xy1.z;
		}

		if (in->data[i].xy2.x>x)
		{
			x=in->data[i].xy2.x;
		}

		if (in->data[i].xy2.y>y)
		{
			y=in->data[i].xy2.y;
		}

		if (in->data[i].xy2.z>z)
		{
			z=in->data[i].xy2.z;
		}
	}

	out->x=x;
	out->y=y;
	out->z=z;

}

void triangles_sub_vec(struct triangles *in,struct vec *v)
{
	int i;
	for (i=0;i<in->len;i++)
	{
		vec_sub(&(in->data[i].xy0),v);
		vec_sub(&(in->data[i].xy1),v);
		vec_sub(&(in->data[i].xy2),v);
	}
}

void triangles_add_vec(struct triangles *in,struct vec *v)
{
	int i;
	for (i=0;i<in->len;i++)
	{
		vec_add(&(in->data[i].xy0),v);
		vec_add(&(in->data[i].xy1),v);
		vec_add(&(in->data[i].xy2),v);
	}
}

void triangles_div_vec(struct triangles *in,struct vec *v)
{
	int i;
	for (i=0;i<in->len;i++)
	{
		vec_div_vec(&(in->data[i].xy0),v);
		vec_div_vec(&(in->data[i].xy1),v);
		vec_div_vec(&(in->data[i].xy2),v);
	}
}

void triangles_mul_vec(struct triangles *in,struct vec *v)
{
	int i;
	for (i=0;i<in->len;i++)
	{
		vec_mul_vec(&(in->data[i].xy0),v);
		vec_mul_vec(&(in->data[i].xy1),v);
		vec_mul_vec(&(in->data[i].xy2),v);
	}
}

void triangles_print(struct triangles *in)
{
	int i;
	for (i=0;i<in->len;i++)
	{
		printf("(%le,%le,%le)\n",in->data[i].xy0.x,in->data[i].xy0.y,in->data[i].xy0.z);
		printf("(%le,%le,%le)\n",in->data[i].xy1.x,in->data[i].xy1.y,in->data[i].xy1.z);
		printf("(%le,%le,%le)\n",in->data[i].xy2.x,in->data[i].xy2.y,in->data[i].xy2.z);
		printf("\n");
	}

	printf("%d triangles\n",in->len);

}



void triangles_add_triangle(struct triangles *obj, double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,int uid,int object_type)
{

	obj->data[obj->len].xy0.x=x0;
	obj->data[obj->len].xy0.y=y0;
	obj->data[obj->len].xy0.z=z0;

	obj->data[obj->len].xy1.x=x1;
	obj->data[obj->len].xy1.y=y1;
	obj->data[obj->len].xy1.z=z1;

	obj->data[obj->len].xy2.x=x2;
	obj->data[obj->len].xy2.y=y2;
	obj->data[obj->len].xy2.z=z2;

	obj->data[obj->len].object_type=object_type;
	obj->data[obj->len].object_uid=uid;
	obj->len++;
	if (obj->len>=obj->max_len)
	{
		obj->max_len+=1000;
		obj->data=(struct triangle *)realloc(obj->data,obj->max_len*sizeof(struct triangle));
		if (obj->data==NULL)
		{
			printf("triangle memory errror\n");
		}
		//printf("%d %d\n",in->triangles,in->triangles_max);
	}

}

void triangles_malloc(struct triangles *tri)
{
	tri->max_len=30;
	tri->data=(struct triangle *)malloc(tri->max_len*sizeof(struct triangle));

	if (tri->data==NULL)
	{
		printf("triangle memory errror\n");
	}
}

void triangles_init(struct triangles *tri)
{
	tri->len=0;
	tri->max_len=0;
	tri->data=NULL;
	tri->edges_calculated=FALSE;

}

void triangles_cal_edges(struct triangles *in)
{
	int i;
	struct triangle *tri;

	for (i=0;i<in->len;i++)
	{
		tri=&(in->data[i]);

		//edge1 = vertex1 - vertex0;
		vec_cpy(&(tri->edge1),&(tri->xy1));
		vec_sub(&(tri->edge1),&(tri->xy0));

		//edge2 = vertex2 - vertex0;
		vec_cpy(&(tri->edge2),&(tri->xy2));
		vec_sub(&(tri->edge2),&(tri->xy0));
	}

	in->edges_calculated=TRUE;

}

double triangles_interpolate(struct triangles *in,struct vec *p)
{
	int i;
	char temp[100];
	struct triangle *tri;
	double ret;
	double max=0.0;
	double m=0.0;
	int i_max=0;
	for (i=0;i<in->len;i++)
	{
		tri=&(in->data[i]);
		if (triangle_vec_within (tri,p)==TRUE)
		{
			m=triangle_get_max_y(tri);
			if (m>max)
			{
				max=m;
				i_max=i;
			}

		}
		//printf("%d %d\n",i,triangle_vec_within (tri,p));
	}

	//sprintf(temp,"tri.dat",i_max);
	//triangle_dump(temp,&(in->data[i_max]));

	//FILE *out=fopen("p.dat","w");
	//fprintf(out,"%le %le %le\n",p->x, p->z,,p->x, p->z));
	//fclose(out);
	ret=triangle_get_y_from_xz(&(in->data[i_max]),p->x, p->z);
	return ret;
}
