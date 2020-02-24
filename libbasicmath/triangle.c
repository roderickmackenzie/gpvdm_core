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


/** @file triangle.c
	@brief Basic low level triangle functions
*/

long double get_sign (struct vec *p0, struct vec *p1, struct vec *p2)
{
    return (p0->x - p2->x) * (p1->z - p2->z) - (p1->x - p2->x) * (p0->z - p2->z);
}

int triangle_vec_within (struct triangle *tri,struct vec *pt)
{
    long double d0, d1, d2;
    int is_neg;
	int is_pos;

    d0 = get_sign(pt, &(tri->xy0), &(tri->xy1));
    d1 = get_sign(pt, &(tri->xy1), &(tri->xy2));
    d2 = get_sign(pt, &(tri->xy2), &(tri->xy0));

    is_neg = (d0 < 0) || (d1 < 0) || (d2 < 0);
    is_pos = (d0 > 0) || (d1 > 0) || (d2 > 0);

    return !(is_neg && is_pos);
}

void triangle_print(struct triangle *in)
{
	printf("(%le,%le,%le)\n",in->xy0.x,in->xy0.y,in->xy0.z);
	printf("(%le,%le,%le)\n",in->xy1.x,in->xy1.y,in->xy1.z);
	printf("(%le,%le,%le)\n",in->xy2.x,in->xy2.y,in->xy2.z);
	printf("\n");

}

void triangle_cog(struct vec *out,struct triangle *in)			//find the center of gravity
{
	out->x=(in->xy0.x+in->xy1.x+in->xy2.x)/3.0;
	out->y=(in->xy0.y+in->xy1.y+in->xy2.y)/3.0;
	out->z=(in->xy0.z+in->xy1.z+in->xy2.z)/3.0;


}


void triangle_norm(struct vec *ret,struct triangle *my_obj)
{
	double mag=0.0;
	struct vec edge0;
	struct vec edge1;
	struct vec n;

	vec_init(&edge0);
	vec_init(&edge1);
	vec_init(&n);


	vec_cpy(&edge0,&(my_obj->xy1));
	vec_sub(&edge0,&(my_obj->xy0));

	vec_cpy(&edge1,&(my_obj->xy2));
	vec_sub(&edge1,&(my_obj->xy0));

	vec_cross(&n,&edge0,&edge1);
	mag=vec_fabs(&n);
	vec_div(&n,mag);

	vec_cpy(ret,&n);
}

double triangle_get_y_from_xz(struct triangle *tri,double x, double z)
{
	double y=0.0;
	struct vec edge0;
	struct vec edge1;
	struct vec n;

	vec_init(&edge0);
	vec_init(&edge1);
	vec_init(&n);

	vec_cpy(&edge0,&(tri->xy1));
	vec_sub(&edge0,&(tri->xy0));

	vec_cpy(&edge1,&(tri->xy2));
	vec_sub(&edge1,&(tri->xy0));

	vec_cross(&n,&edge0,&edge1);

	y=(n.x*(x-tri->xy0.x)+n.z*(z-tri->xy0.z))/(-n.y)+tri->xy0.y;
	return y;
}

void triangle_dump(char *file_name,struct triangle *tri)
{
	FILE *out;
	out=fopen(file_name,"w");
	fprintf(out,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
	fprintf(out,"%le %le %le\n",tri->xy1.z,tri->xy1.x,tri->xy1.y);
	fprintf(out,"%le %le %le\n",tri->xy2.z,tri->xy2.x,tri->xy2.y);
	fprintf(out,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
	fclose(out);

}


double triangle_get_min_y(struct triangle* tri)
{
	double min=tri->xy0.y;
	if (min>tri->xy1.y)
	{
		min=tri->xy1.y;
	}

	if (min>tri->xy2.y)
	{
		min=tri->xy2.y;
	}

	return min;
}

double triangle_get_max_y(struct triangle* tri)
{
	double min=tri->xy0.y;
	if (min<tri->xy1.y)
	{
		min=tri->xy1.y;
	}

	if (min<tri->xy2.y)
	{
		min=tri->xy2.y;
	}

	return min;
}

void ray_tri_flip_y_axis(struct triangle* tri,double y)
{
	vec_flip_y_axis(&(tri->xy0),y);
	vec_flip_y_axis(&(tri->xy1),y);
	vec_flip_y_axis(&(tri->xy2),y);
}

void ray_tri_sub_y(struct triangle* tri,double y)
{
	struct vec a;
	vec_set(&a,0,y,0);
	vec_sub(&(tri->xy0),&a);
	vec_sub(&(tri->xy1),&a);
	vec_sub(&(tri->xy2),&a);

}

void ray_tri_add_y(struct triangle* tri,double y)
{
	struct vec a;
	vec_set(&a,0,y,0);
	vec_add(&(tri->xy0),&a);
	vec_add(&(tri->xy1),&a);
	vec_add(&(tri->xy2),&a);

}
