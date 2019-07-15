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

/** @file ray.c
	@brief Ray tracing for the optical model, this should really be split out into it's own library.
*/

void image_init(struct image *in)
{
	in->lines=0;
	in->nrays=0;
	in->objects=0;
	in->n_start_rays=0;
	in->ray_wavelength_points=0;
	in->extract_eff=NULL;
	in->ray_auto_run=FALSE;
	in->escape_bins=0;
	in->ray_xsrc=-1.0;
	in->ray_ysrc=-1.0;
	in->ray_zsrc=-1.0;
	in->ray_theta_start=0.0;
	in->ray_theta_stop=360.0;
}

int between(double v, double x0, double x1)
{
	double min=0.0;
	double max=0.0;
	if (x0>x1)
	{
		min=x1;
		max=x0;
	}else
	{
		min=x0;
		max=x1;
	}

	if ((v-min)>=-1e-12)
	{
		if ((v-max)<=1e-12)
		{
			return 0;
		}
	}

	return -1;
}


void ray_reset(struct image *in)
{
	in->nrays=0;
}

void add_ray(struct simulation *sim,struct image *in,struct vec *start,struct vec *dir,double mag)
{
	if (mag>1e-3)
	{
		vec_cpy(&(in->rays[in->nrays].xy),start);
		vec_cpy(&(in->rays[in->nrays].dir),dir);
		in->rays[in->nrays].state=WAIT;
		in->rays[in->nrays].bounce=0;
		in->rays[in->nrays].mag=mag;
		in->nrays++;
		if (in->nrays>=RAY_MAX)
		{
			printf_log(sim,"too many rays!\n");
			exit(0);
		}
	}
		
}




void obj_norm(struct vec *ret,struct triangle *my_obj)
{
	double norm=0.0;
	vec_cpy(ret,&(my_obj->xy0));
	vec_sub(ret,&(my_obj->xy1));
	vec_swap(ret);
	norm=vec_fabs(ret);
	vec_div(ret,norm);
}


int ray_intersect(struct vec *ret,struct triangle *my_obj,struct ray *my_ray)
{
double x1=my_obj->xy0.x;
double x2=my_obj->xy1.x;
double x3=my_ray->xy.x;
double x4=my_ray->xy.x+my_ray->dir.x;

double y1=my_obj->xy0.y;
double y2=my_obj->xy1.y;
double y3=my_ray->xy.y;
double y4=my_ray->xy.y+my_ray->dir.y;

double x12 = x1 - x2;

double x34 = x3 - x4;
double y12 = y1 - y2;
double y34 = y3 - y4;

double c = x12 * y34 - y12 * x34;
double x=0.0;
double y=0.0;
int hit=0;
if (fabs(c) < 1e-10)
{
  
}
else
{
	double a = x1 * y2 - y1 * x2;
	double b = x3 * y4 - y3 * x4;

	x = (a * x34 - b * x12) / c;
	y = (a * y34 - b * y12) / c;

	//printf_log(sim,"here  x0=%le x=%le %le x1=%le %d\n",my_obj->xy0.x,x,my_obj->xy1.x,my_ray->xy.x,between(x,my_obj->xy0.x,my_obj->xy1.x));
	//getchar();
	if (between(x,my_obj->xy0.x,my_obj->xy1.x)==0)
	{

		if (between(y,my_obj->xy0.y,my_obj->xy1.y)==0)
		{

			if (my_ray->dir.y>=0)
			{
				if (y>my_ray->xy.y)
				{
					hit=1;
				}else
				{
					hit=0;
				}
			}

			if (my_ray->dir.y<=0)
			{
				if (y<my_ray->xy.y)
				{
					hit=1;
				}else
				{
					hit=0;
				}
			}


			if (my_ray->dir.x>0)
			{
				if (x>my_ray->xy.x)
				{
					hit=1;
				}else
				{
					hit=0;
				}
			}

			if (my_ray->dir.x<0)
			{
				if (x<my_ray->xy.x)
				{
					hit=1;
				}else
				{
					hit=0;
				}
			}
		}
	}
}


ret->x=x;
ret->y=y;

return hit;
}

int search_ojb(struct image *in,int ray)
{
int i;
int found=-1;
int pos=-1;

double min_dist=1000.0;
double dist=0;

struct vec ret;
vec_init(&ret);

struct vec tmp;
vec_init(&tmp);

struct vec store;
vec_init(&store);

	for (i=0;i<in->lines;i++)
	{
		found=ray_intersect(&ret,&(in->p[i]),&(in->rays[ray]));

		if (found==1)
		{
			vec_cpy(&tmp,&ret);
			//vec_print(&ret);
			//getchar();
			//vec_print(&(in->rays[ray].xy));			
			vec_sub(&tmp,&(in->rays[ray].xy));
			dist=vec_fabs(&tmp);

			if (dist>1e-12)
			{
				if (dist<min_dist)
				{
					pos=i;
					vec_cpy(&(in->rays[ray].xy_end),&ret);
					min_dist=dist;
				}
			}
		}
		
	}


return pos;	
}

int activate_rays(struct image *in)
{
	int i=0;
	int changed=0;
	for (i=0;i<in->nrays;i++)
	{
		if (in->rays[i].state==WAIT)
		{
			in->rays[i].state=READY;
			changed++;
		}
		
	}

return changed;
}

int pnpoly(struct image *in, struct vec *xy,int id)
{
	int c=1;
	int i=0;
	double x=xy->x;
	double y=xy->y;
	for (i=0;i<in->lines;i++)
	{
		if (in->p[i].id==id)
		{
			if (between(y,in->p[i].xy0.y,in->p[i].xy1.y)==0)
			{
				if (x < (in->p[i].xy1.x - in->p[i].xy0.x) * (y - in->p[i].xy0.y) / (in->p[i].xy1.y - in->p[i].xy0.y) + in->p[i].xy0.x)
				{
					c = !c;
				}
					
			}
		}
		
	}
	
	return c;
}

void get_refractive(struct image *in,double *alpha,double *n0,double *n1,int ray)
{
	struct vec tmp;
	vec_init(&tmp);

	struct vec back;
	vec_init(&back);

	struct vec fwd;
	vec_init(&fwd);

	vec_cpy(&tmp,&(in->rays[ray].dir));
	vec_mul(&tmp,-1e-10);
	vec_cpy(&back,&(in->rays[ray].xy_end));
	vec_add(&back,&tmp);

	vec_cpy(&tmp,&(in->rays[ray].dir));
	vec_mul(&tmp,1e-10);
	vec_cpy(&fwd,&(in->rays[ray].xy_end));
	vec_add(&fwd,&tmp);
	//vec_print(&fwd);
	//vec_print(&back);

	int i;
	int i_fwd=-1;
	int i_back=-1;
	
	for (i=in->objects-1;i>=0;i--)
	{
		if (pnpoly(in, &back,i)==0)
		{
			i_back=i;
			break;
		}
	}

	for (i=in->objects-1;i>=0;i--)
	{
		if (pnpoly(in, &fwd,i)==0)
		{
			i_fwd=i;
			break;
		}
	}
	

	if ((i_fwd!=-1))
	{
		*n1=in->obj_n[i_fwd];
	}

	if (i_back!=-1)
	{
		*n0=in->obj_n[i_back];
		*alpha=in->obj_alpha[i_back];
	}

}
    
int propergate_next_ray(struct simulation *sim,struct image *in)
{
	struct vec n;
	vec_init(&n);

	struct vec n_inv;
	vec_init(&n_inv);
	
	struct vec r;
	vec_init(&r);

	struct vec t;
	vec_init(&t);

	struct vec temp;
	vec_init(&temp);

	double threshold=0.0;
	double ang_out=0.0;
	int ray=0;

	double R=0.0;
	double T=0.0;
	double mag=0.0;
	
	for (ray=0;ray<in->nrays;ray++)
	{
		if (in->rays[ray].state==READY)
		{
			
			in->rays[ray].state=DONE;

			int item=search_ojb(in,ray);

			if (item==-1)
			{
				vec_cpy(&in->rays[ray].xy_end,&in->rays[ray].xy);
				
			}else
			{
				double dist=vec_dist(&(in->rays[ray].xy),&(in->rays[ray].xy_end));

				int bounce=in->rays[ray].bounce;
				mag=in->rays[ray].mag;
				bounce=bounce+1;

				double n0=1.0;
				double n1=1.0;
				double alpha=1e9;
				get_refractive(in,&alpha,&n0,&n1,ray);

				//Calculate norm of object
				obj_norm(&n,&(in->p[item]));
				vec_cpy(&n_inv,&n);
				vec_mul(&n_inv,-1.0);

				//vec_print(&n);

				////Cal normal to surface
				double a=vec_overlap(&n,&(in->rays[ray].dir));
				vec_mul(&n,-1.0);
				vec_mul(&n_inv,-1.0);
				double b=vec_overlap(&n,&(in->rays[ray].dir));

				if (a<b)
				{
					vec_mul(&n,-1.0);
					vec_mul(&n_inv,-1.0);

				}

				/////////////Calculate reflected ray.
				double ret=2.0*vec_dot(&(in->rays[ray].dir),&n);

				vec_cpy(&temp,&n);
				vec_mul(&temp,ret);

				vec_cpy(&r,&(in->rays[ray].dir));
				vec_sub(&r,&temp);
				//////Angle between inceident ray and surface
				double dot=0.0;
				double ang_in=0.0;
				dot=vec_dot(&n_inv,&(in->rays[ray].dir));
				ang_in=PI/2.0-acos(dot);



				if (n1>=n0)
				{
					threshold=0.0;
				}else
				{
					threshold=asin(n1/n0);
				}

				///////////////Calculate Snell's law in vector form.
				struct vec left;
				vec_init(&left);
				vec_cross(&left,&n,&(in->rays[ray].dir));
				vec_mul(&left,-1.0);
				vec_cross(&left,&n,&left);
				vec_mul(&left,n0/n1);


				struct vec cr;
				vec_init(&cr);
				vec_cross(&cr,&n,&(in->rays[ray].dir));

				dot=0.0;
				dot=vec_dot(&cr,&cr);
				dot=dot*(pow((n0/n1),2.0));
				double radicand=1.0-dot;

				if (radicand>=0.0)
				{
					dot=sqrt(radicand);

					struct vec right;
					vec_init(&right);
					vec_cpy(&right,&n);
					vec_mul(&right,dot);
					vec_sub(&left,&right);
					
					vec_cpy(&t,&left);
				}else
				{
					vec_set(&t,0.0,0.0,0.0);
				}

				//////Out angle

				if ((ang_in>threshold)&&(radicand>=0.0))
				{
					dot=vec_dot(&n_inv,&t);
					ang_out=PI/2.0-acos(dot);

					R=((n0*cos(PI/2.0-ang_in)-n1*cos(PI/2.0-ang_out))/(n0*cos(PI/2.0-ang_in)+n1*cos(PI/2.0-ang_out)));
					R=fabs(R);
					T=1.0-R;

				}else
				{
					ang_out=-100.0;
					R=1.0;
					T=0.0;
				}

				//////////////////////////////


				//vec_print(&r);

				if (bounce<100)
				{
					if (in->p[item].edge==FALSE)
					{
						double abs=exp(-alpha*dist);

						add_ray(sim,in,&(in->rays[ray].xy_end),&r,R*mag*abs);
						in->rays[in->nrays-1].bounce=bounce;


						if (ang_in>=threshold)
						{
							add_ray(sim,in,&(in->rays[ray].xy_end),&(t),T*mag);
							in->rays[in->nrays-1].bounce=bounce;
						}	
							
					}
				}	
			}
		}
	}

return 0;
}


