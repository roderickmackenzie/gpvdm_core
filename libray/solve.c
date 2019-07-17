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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dump_ctrl.h>

/** @file solve.c
	@brief This will call the ray tracer for the standard case.
*/

double get_rand()
{
	double r=0.0;
	r = rand();
	r=(double)r/(double)RAND_MAX;
	return r;
}

void ray_solve(struct simulation *sim,struct device *in, int l)
{
	ray_escape_angle_reset(&(in->my_image),l);

	int i;
	double x_vec;
	double y_vec;
	double z_vec;

	long double lam =in->my_image.lam[l];
	int x=0;
	double theta=0.0;
	double dtheta=0.0;

	double phi=0.0;
	double dphi=0.0;

	double eff=0.0;
	int sims=0;
	long double mag=0.0;

	char name[400];
	char out_dir[400];
	char one[200];
	struct stat st = {0};

	FILE *out;
	sprintf(out_dir,"%s/ray_trace",get_output_path(sim));
	sprintf(name,"%s/light_ray_%d.dat",out_dir,(int)(lam*1e9));

	if (stat(out_dir, &st) == -1)
	{
		mkdir(out_dir, 0700);
	}
	
	if ((out=fopen(name,"w")))
	{
		fclose(out);
	}
	
	//if (get_dump_status(sim,dump_optics)==TRUE)
	//{

	mag=inter_get_hard(&(in->my_image.input_spectrum),lam);

	sprintf(one,"Ray tracing at %Lf nm mag=%Lf\n",lam*1e9,mag);
	//printf("%Lf\n",in->my_image.input_spectrum.data[in->my_image.input_spectrum.len-1]);
	//inter_dump(sim,&(in->my_image.input_spectrum));
	waveprint(sim,one,lam*1e9);
	//}
	light_update_ray_mat(sim,&(in->my_epitaxy),&(in->my_image),lam);

	//for (x=0;x<in->my_image.n_start_rays;x++)


	x=0;
	phi=in->my_image.ray_phi_start;
	dphi=((in->my_image.ray_phi_stop-in->my_image.ray_phi_start)/(double)in->my_image.phi_steps);

	while(phi<in->my_image.ray_phi_stop)
	{
		theta=in->my_image.ray_theta_start;
		dtheta=((in->my_image.ray_theta_stop-in->my_image.ray_theta_start)/(double)in->my_image.theta_steps);

		while(theta<in->my_image.ray_theta_stop)
		{
			//theta=70.0;//+10.0;
			//theta=get_rand()*360.0;
			//printf("%lf %lf\n",dtheta,theta);
			//getchar();
			x_vec=cos(2*PI*(theta/360.0));
			y_vec=sin(2*PI*(theta/360.0));
			z_vec=sin(2*PI*(phi/360.0));

			struct vec start;
			vec_init(&start);

			struct vec dir;
			vec_init(&dir);

			vec_set(&start,in->my_image.start_rays[x].x,in->my_image.start_rays[x].y,in->my_image.start_rays[x].z);
			vec_set(&dir,x_vec,y_vec,z_vec);
			
			add_ray(sim,&in->my_image,&start,&dir,mag);
			activate_rays(&in->my_image);

			//dump_plane(sim,&in->my_image);
			//dump_plane_to_file(name,&in->my_image);

			//vec_print(&start);
			//getchar();
			//exit(0);
			int ret=0;
			for (i=0;i<50;i++)
			{
				propergate_next_ray(sim,&in->my_image);
				//dump_plane(sim,&in->my_image);
				//dump_plane_to_file(name,&in->my_image);
				//getchar();
				ret=activate_rays(&in->my_image);

				//dump_plane(sim,&in->my_image);
				//printf("bounce\n");
				//getchar();

				if (ret==0)
				{
					break;
				}

			}
			//getchar();
			//exit(0);
			//if (get_dump_status(sim,dump_ray_trace_map)==TRUE)
			//{
			//getchar();
			//}
			double tot=ray_cal_escape_angle(&(in->my_image),l);
			//printf("%lf %lf\n",theta,tot);

			double e=get_eff(&in->my_image);
			eff+=e;
			sims++;
			//ray_reset(&in->my_image);

			theta+=dtheta;
			//getchar();
			//exit(0);
		}

		phi+=dphi;
		printf("%lf\n",phi);
		//getchar();
	}

	dump_plane_to_file(sim,"lines.dat",&in->my_image);

	in->my_image.extract_eff[l]=eff/((double)sims);

}

