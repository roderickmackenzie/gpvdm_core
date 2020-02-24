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
#include <string.h>
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
#include <triangles.h>

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

	strcpy(li->light_profile,"box");
	triangles_init((&(li->light_profile_tri)));

	dim_light_init(dim);
	matrix_init(&(li->mx));
}

