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


/** @file light_memory.c
	@brief Deals with memory allocation for the light model.
*/

#include <stdlib.h>
#include "util.h"
#include "gpvdm_const.h"
#include "light.h"
#include "device.h"
#include "config.h"
#include "util.h"
#include "lang.h"
#include "log.h"
#include <light_fun.h>
#include <matrix.h>
#include <memory.h>
#include <device_fun.h>


static int unused __attribute__((unused));

void light_init(struct light *li)
{
	struct dim_light *dim=&(li->dim);
	li->last_Psun= -1000.0;
	li->last_laser_eff= -1000.0;
	li->last_wavelength_laser= -1000.0;
	li->laser_wavelength= -1.0;
	li->laser_pos= -1;
	li->laser_wavelength= -1.0;
	li->lstart= -1.0;
	li->spotx= -1.0;
	li->spoty= -1.0;
	li->pulseJ= -1.0;
	li->pulse_width= -1.0;
	li->G_percent=NULL;
	li->disable_cal_photon_density=FALSE;
	li->finished_solveing=FALSE;
	li->print_wavlengths=TRUE;
	li->save_data_to_disk=TRUE;

	li->Ep=NULL;
	li->Epz=NULL;
	li->En=NULL;
	li->Enz=NULL;
	li->n=NULL;
	li->alpha0=NULL;
	li->alpha=NULL;
	li->photons=NULL;
	li->photons_asb=NULL;
	li->pointing_vector=NULL;
	li->E_tot_r=NULL;
	li->E_tot_i=NULL;
	li->H=NULL;
	li->Gn=NULL;
	li->Gp=NULL;

	//long double complex
	li->t=NULL;
	li->r=NULL;
	li->nbar=NULL;

	li->Htot=NULL;
	li->photons_tot=NULL;

	li->sun_E=NULL;
	li->sun=NULL;
	li->sun_norm=NULL;
	li->sun_photons=NULL;
	li->reflect=NULL;
	li->transmit=NULL;
	li->extract_eff=NULL;

	dim_light_init(dim);
	matrix_init(&(li->mx));
}

void light_malloc(struct simulation *sim,struct light *li)
{
printf_log(sim,"alloc: light_memory\n");
	int l=0;
	struct dim_light *dim=&(li->dim);
	struct matrix *mx=&(li->mx);

	dim_light_malloc(dim);

	//long double zxyl
	malloc_light_zxyl_long_double(dim,&(li->Ep));
	malloc_light_zxyl_long_double(dim,&(li->Epz));
	malloc_light_zxyl_long_double(dim,&(li->En));
	malloc_light_zxyl_long_double(dim,&(li->Enz));
	malloc_light_zxyl_long_double(dim,&(li->n));
	malloc_light_zxyl_long_double(dim,&(li->alpha0));
	malloc_light_zxyl_long_double(dim,&(li->alpha));
	malloc_light_zxyl_long_double(dim,&(li->photons));
	malloc_light_zxyl_long_double(dim,&(li->photons_asb));
	malloc_light_zxyl_long_double(dim,&(li->pointing_vector));
	malloc_light_zxyl_long_double(dim,&(li->E_tot_r));
	malloc_light_zxyl_long_double(dim,&(li->E_tot_i));
	malloc_light_zxyl_long_double(dim,&(li->H));


	//long double zxy
	malloc_light_zxy_long_double(dim,&(li->Gn));
	malloc_light_zxy_long_double(dim,&(li->Gp));
	malloc_light_zxy_long_double(dim,&(li->Htot));
	malloc_light_zxy_long_double(dim,&(li->photons_tot));

	//long double complex
	malloc_light_zxyl_long_double_complex(dim,&(li->t));
	malloc_light_zxyl_long_double_complex(dim,&(li->r));
	malloc_light_zxyl_long_double_complex(dim,&(li->nbar));

	//zxy_p_object
	malloc_light_zxy_p_object(dim, &(li->obj));

	malloc_light_l_long_double(dim,&(li->sun_E));
	malloc_light_l_long_double(dim,&(li->sun));
	malloc_light_l_long_double(dim,&(li->sun_norm));
	malloc_light_l_long_double(dim,&(li->sun_photons));
	malloc_light_l_long_double(dim,&(li->reflect));
	malloc_light_l_long_double(dim,&(li->transmit));
	malloc_light_l_long_double(dim,&(li->extract_eff));

	for (l=0;l<dim->llen;l++)
	{
		li->extract_eff[l]=1.0;
	}

	//sort out the matrix
	mx->nz=0.0;
	mx->nz+=dim->ylen;
	mx->nz+=dim->ylen-1;
	mx->nz+=dim->ylen;

	mx->nz+=dim->ylen-1;
	mx->nz+=dim->ylen; //t
	mx->nz+=dim->ylen-1;
	//li->N+=1;
	mx->M=dim->ylen+dim->ylen;
	mx->complex_matrix=TRUE;
	matrix_malloc(sim,mx);

}



void light_free_memory(struct simulation *sim,struct light *li)
{
	struct dim_light *dim=&(li->dim);
	struct matrix *mx=&(li->mx);

	//long double zxyl
	free_light_zxyl_long_double(dim,&(li->Ep));
	free_light_zxyl_long_double(dim,&(li->Epz));
	free_light_zxyl_long_double(dim,&(li->En));
	free_light_zxyl_long_double(dim,&(li->Enz));
	free_light_zxyl_long_double(dim,&(li->n));
	free_light_zxyl_long_double(dim,&(li->alpha0));
	free_light_zxyl_long_double(dim,&(li->alpha));
	free_light_zxyl_long_double(dim,&(li->photons));
	free_light_zxyl_long_double(dim,&(li->photons_asb));
	free_light_zxyl_long_double(dim,&(li->pointing_vector));
	free_light_zxyl_long_double(dim,&(li->E_tot_r));
	free_light_zxyl_long_double(dim,&(li->E_tot_i));
	free_light_zxyl_long_double(dim,&(li->H));

	//long double zxy
	free_light_zxy_long_double(dim,&(li->Gn));
	free_light_zxy_long_double(dim,&(li->Gp));
	free_light_zxy_long_double(dim,&(li->Htot));
	free_light_zxy_long_double(dim,&(li->photons_tot));

	//long double complex
	free_light_zxyl_long_double_complex(dim,&(li->t));
	free_light_zxyl_long_double_complex(dim,&(li->r));
	free_light_zxyl_long_double_complex(dim,&(li->nbar));

	//zxy_p_object
	free_light_zxy_p_object(dim, &(li->obj));

	//long double
	free_light_l_long_double(dim,&(li->sun_E));
	free_light_l_long_double(dim,&(li->sun));
	free_light_l_long_double(dim,&(li->sun_norm));
	free_light_l_long_double(dim,&(li->sun_photons));
	free_light_l_long_double(dim,&(li->reflect));
	free_light_l_long_double(dim,&(li->transmit));
	free_light_l_long_double(dim,&(li->extract_eff));

/////////////
	matrix_free(sim,mx);


	inter_free(&(li->sun_read));
	dim_light_free(dim);

	printf_log(sim,_("Freeing memory from the optical model\n"));
}
