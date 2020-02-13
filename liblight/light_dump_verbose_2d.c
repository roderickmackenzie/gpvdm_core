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

/** @file light_dump_vervose_2d.c
	@brief Dumps the optical fields in 2D.
*/


#include <string.h>
#include <sys/stat.h>
#include "util.h"
#include "gpvdm_const.h"
#include "dump_ctrl.h"
#include "light.h"
#include "dat_file.h"
#include <cal_path.h>
#include <lang.h>

void light_dump_verbose_2d(struct simulation *sim,struct light *li)
{
	FILE *out;
	int x=0;
	int y=0;
	int z=0;
	int l=0;
	struct dat_file buf;
	char line[1024];
	char temp[1024];
	struct dim_light *dim=&li->dim;
	struct epitaxy *epi=li->epi;
	long double device_start =epi->device_start;
	buffer_init(&buf);

	out=fopena(li->dump_dir,"light_2d_Ep.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			fprintf(out,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,gpow(gpow(li->Ep[z][x][y][l],2.0)+gpow(li->Epz[z][x][y][l],2.0),0.5));

		}

	fprintf(out,"\n");
	}
	fclose(out);

	out=fopena(li->dump_dir,"light_2d_En.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			fprintf(out,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,gpow(gpow(li->En[z][x][y][l],2.0)+gpow(li->Enz[z][x][y][l],2.0),0.5));
		}

	fprintf(out,"\n");
	}
	fclose(out);

	out=fopena(li->dump_dir,"light_2d_E_mod.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			fprintf(out,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,gpow(gpow(li->Ep[z][x][y][l]+li->En[z][x][y][l],2.0)+gpow(li->Enz[z][x][y][l]+li->Epz[z][x][y][l],2.0),1.0));
		}

	fprintf(out,"\n");
	}
	fclose(out);



	buffer_malloc(&buf);
	buf.y_mul=1e9;
	buf.x_mul=1e9;
	strcpy(buf.title,"Refractive index (real) coefficient");
	strcpy(buf.type,"heat");
	strcpy(buf.x_label,"Position");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Density");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"m^{-3}");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=dim->llen;
	buf.y=dim->ylen;
	buf.z=1;
	buffer_add_info(sim,&buf);

	sprintf(temp,"#data\n");
	buffer_add_string(&buf,temp);

	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(line,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,li->n[z][x][y][l]);
			buffer_add_string(&buf,line);
		}

	buffer_add_string(&buf,"\n");
	}

	sprintf(temp,"#end\n");
	buffer_add_string(&buf,temp);

	buffer_dump_path(sim,li->dump_dir,"2d_n.dat",&buf);
	buffer_free(&buf);


	out=fopena(li->dump_dir,"light_lambda_sun.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		fprintf(out,"%Le %Le\n",dim->l[l],li->sun[l]);
	}
	fclose(out);

	out=fopena(li->dump_dir,"light_lambda_sun_norm.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		fprintf(out,"%Le %Le\n",dim->l[l],li->sun_norm[l]);
	}
	fclose(out);

	out=fopena(li->dump_dir,"light_lambda_sun_photons.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		fprintf(out,"%Le %Le\n",dim->l[l],li->sun_photons[l]);
	}
	fclose(out);

	buffer_malloc(&buf);
	buf.y_mul=1e9;
	buf.x_mul=1e9;
	strcpy(buf.title,"Optical absorption coefficient");
	strcpy(buf.type,"heat");
	strcpy(buf.x_label,"Position");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Density");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"m^{-3}");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=dim->llen;
	buf.y=dim->ylen;
	buf.z=1;
	buffer_add_info(sim,&buf);

	sprintf(temp,"#data\n");
	buffer_add_string(&buf,temp);

	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(line,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,li->alpha[z][x][y][l]);
			buffer_add_string(&buf,line);
		}

	buffer_add_string(&buf,"\n");
	}

	sprintf(temp,"#end\n");
	buffer_add_string(&buf,temp);

	buffer_dump_path(sim,li->dump_dir,"2d_alpha.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1e9;
	buf.x_mul=1e9;
	strcpy(buf.title,"Optical absorption coefficient");
	strcpy(buf.type,"heat");
	strcpy(buf.x_label,"Position");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Density");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"m^{-3}");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=dim->llen;
	buf.y=dim->ylen;
	buf.z=1;
	buffer_add_info(sim,&buf);

	sprintf(temp,"#data\n");
	buffer_add_string(&buf,temp);

	for (l=0;l<dim->llen;l++)
	{
		for (y=0;y<dim->ylen;y++)
		{
			sprintf(line,"%Le %Le %Le\n",dim->l[l],dim->y[y]-device_start,li->n[z][x][y][l]);
			buffer_add_string(&buf,line);
		}

	buffer_add_string(&buf,"\n");
	}

	sprintf(temp,"#end\n");
	buffer_add_string(&buf,temp);

	buffer_dump_path(sim,li->dump_dir,"light_lambda_n.dat",&buf);
	buffer_free(&buf);

	out=fopena(li->dump_dir,"light_sun_wavelength_E.dat","w");
	for (l=0;l<dim->llen;l++)
	{
		fprintf(out,"%Le %Le\n",dim->l[l],li->sun_E[l]);
	}
	fclose(out);


}
