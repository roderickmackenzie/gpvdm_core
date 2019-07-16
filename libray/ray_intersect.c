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
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <ray_fun.h>

/** @file ray_intersect.c
	@brief Ray intersect
*/

//bool RayIntersectsTriangle(Vector3D rayOrigin, Vector3D rayVector, Triangle* inTriangle, Vector3D& outIntersectionPoint)
int ray_intersect(struct vec *ret,struct triangle *my_obj,struct ray *my_ray)
{
    double EPSILON = 1e-18;
	double a,f,u,v,t;
    //Vector3D vertex0 = inTriangle->vertex0;
    //Vector3D vertex1 = inTriangle->vertex1;  
    //Vector3D vertex2 = inTriangle->vertex2;
    //Vector3D edge1, edge2, h, s, q;
	struct vec h;
	struct vec s;
	struct vec q;
	struct vec edge1;
	struct vec edge2;
	vec_init(&edge1);
	vec_init(&edge2);
	vec_init(&h);
	vec_init(&s);
	vec_init(&q);



	//edge1 = vertex1 - vertex0;
	vec_cpy(&edge1,&(my_obj->xy1));
	vec_sub(&edge1,&(my_obj->xy0));

	//edge2 = vertex2 - vertex0;
	vec_cpy(&edge2,&(my_obj->xy2));
	vec_sub(&edge2,&(my_obj->xy0));

	//h = rayVector.crossProduct(edge2);
	//vec_print(&(edge1));
	//vec_print(&(edge2));

	vec_cross(&h,&(my_ray->dir),&edge2);

	//vec_print(&(h));

    //a = edge1.dotProduct(h);
	a=vec_dot(&edge1,&h);
	//printf("%e\n",a);
    if (a > -EPSILON && a < EPSILON)
	{
		//printf("exit1\n");
		return FALSE;    // This ray is parallel to this triangle.
	}
    f = 1.0/a;

	//s = rayOrigin - vertex0;
	vec_cpy(&s,&(my_ray->xy));
	vec_sub(&s,&(my_obj->xy0));

    //u = f * s.dotProduct(h);
	u=f*vec_dot(&s,&h);
  
	if (u < 0.0 || u > 1.0)
	{
		//printf("exit2\n");
		return FALSE;
	}

    //q = s.crossProduct(edge1);
	vec_cross(&q,&(s),&edge1);

	//v = f * rayVector.dotProduct(q);
	v=f*vec_dot(&(my_ray->dir),&q);

    if (v < 0.0 || (u + v) > 1.0)
	{
		//printf("exit3\n");
		return FALSE;
	}

    // At this stage we can compute t to find out where the intersection point is on the line.
    //float t = f * edge2.dotProduct(q);
	t= f * vec_dot(&edge2,&q);
    if (t > EPSILON) // ray intersection
    {
		//outIntersectionPoint = rayOrigin + rayVector * t;
		vec_cpy(ret,&(my_ray->dir));
		vec_mul(ret,t);
		vec_add(ret,&(my_ray->xy));
        return TRUE;
    }else // This means that there is a line intersection but not a ray intersection.
    {
		//printf("exit4\n");
		return FALSE;
	}
}

