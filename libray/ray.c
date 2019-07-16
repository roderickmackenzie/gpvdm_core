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
	if (mag>5e-9)
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




void triangle_norm(struct vec *ret,struct triangle *my_obj)
{
	double norm=0.0;
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


void get_refractive(struct simulation *sim,struct image *in,double *alpha,double *n0,double *n1,struct ray *my_ray)
{
	struct vec tmp;
	vec_init(&tmp);

	struct ray back;
	vec_cpy(&(back.dir),&(my_ray->dir));
	vec_cpy(&tmp,&(my_ray->dir));
	vec_mul(&tmp,1e-10);
	vec_cpy(&(back.xy),&(my_ray->xy_end));
	vec_sub(&(back.xy),&(tmp));

	struct ray fwd;
	vec_cpy(&(fwd.dir),&(my_ray->dir));
	vec_cpy(&tmp,&(my_ray->dir));
	vec_mul(&tmp,1e-10);
	vec_cpy(&(fwd.xy),&(my_ray->xy_end));
	vec_add(&(fwd.xy),&(tmp));


	int i;
	int i_fwd=-1;
	int i_back=-1;
	int obj=0;
	i_back=search_object(sim,in,&back);	
	i_fwd=search_object(sim,in,&fwd);
	dump_plane_to_file("lines.dat",in);

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

			int item=search_triangle(sim,in,&in->rays[ray]);
			
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
				
				get_refractive(sim,in,&alpha,&n0,&n1,&(in->rays[ray]));

				//Calculate norm of the surface
				triangle_norm(&n,&(in->p[item]));

				vec_cpy(&n_inv,&n);
				vec_mul(&n_inv,-1.0);

				//The normal will not be relative to the incident ray therefore
				//change the sign until it points in the same direction as the incident ray
				double a=vec_overlap(&n,&(in->rays[ray].dir));
				vec_mul(&n,-1.0);
				vec_mul(&n_inv,-1.0);
				double b=vec_overlap(&n,&(in->rays[ray].dir));

				if (a>b)
				{
					vec_mul(&n,-1.0);
					vec_mul(&n_inv,-1.0);

				}

				//vec_print(&n);
				//vec_print(&n_inv);

				/////////////Calculate reflected ray.
				//https://math.stackexchange.com/questions/2235997/reflecting-ray-on-triangle-in-3d-space
				double ret=2.0*vec_dot(&(in->rays[ray].dir),&n);

				vec_cpy(&temp,&n);
				vec_mul(&temp,ret);

				vec_cpy(&r,&(in->rays[ray].dir));
				vec_sub(&r,&temp);
				
				//////Angle between incident ray and surface
				double dot=0.0;
				double ang_in=0.0;
				dot=vec_dot(&n_inv,&(in->rays[ray].dir));
				ang_in=vec_ang(&n,&(in->rays[ray].dir));

				
				if (n1>=n0)
				{
					threshold=0.0;
				}else
				{
					threshold=asin(n1/n0);
				}

				printf("threshold=%lf\n",360.0*threshold/(2.0*3.1415));
				printf("ang in=%lf\n",360.0*ang_in/(2.0*3.1415));
				printf("%lf\n",n0);
				printf("%lf\n",n1);

				//Vector form of snell's law taken from Wikipedia
				R=0.5;
				T=0.5;
				double sr=n0/n1;
				double sc=-1.0*vec_dot(&n,&(in->rays[ray].dir));
				vec_cpy(&t,&(in->rays[ray].dir));
				vec_mul(&t,sr);
				double inner=sr*sc-sqrt(1.0-sr*sr*(1.0-sc*sc));
				vec_cpy(&temp,&n);
				vec_mul(&temp,inner);
				vec_sub(&t,&temp);

				//If we are not fully internally reflecting then we need to calculate the tx/rx proportions
				
				if (ang_in<threshold)
				{
					ang_out=vec_ang(&n,&t);
					R=((n0*cos(ang_in)-n1*cos(ang_out))/(n0*cos(ang_in)+n1*cos(ang_out)));
					R=fabs(R);
					T=1.0-R;

				}else
				{
					ang_out=-100.0;
					R=1.0;
					T=0.0;
				}

				//printf("angle in=%lf\nangle out=%lf\n",360.0*ang_in/(2.0*3.1415),360.0*ang_out/(2.0*3.1415));
				//printf("T=%lf\nR=%lf %lf\n",T,R,((n0*cos(0.0)-n1*cos(0.0))/(n0*cos(0.0)+n1*cos(0.0))));

				//getchar();
				//////////////////////////////


				//vec_print(&r);

				if (bounce<100)
				{
					printf("%d\n",in->p[item].edge);
					if (in->p[item].edge==FALSE)
					{
						double abs=exp(-alpha*dist);
					
						add_ray(sim,in,&(in->rays[ray].xy_end),&r,R*mag*abs);
						in->rays[in->nrays-1].bounce=bounce;
						
						if (ang_in<threshold)
						{
							//printf("prop>>>>");
							//vec_print(&t);
							//getchar();
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


