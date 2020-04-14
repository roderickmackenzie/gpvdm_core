//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 


/** @file device.c
	@brief Initialize the device structure.
*/

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <device.h>
#include <string.h>
#include <dump.h>
#include <mesh.h>
#include <ray_fun.h>
#include <newton_tricks.h>
#include <memory.h>
#include <circuit.h>
#include <shape.h>
#include <heat.h>
#include <heat_fun.h>


static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));

void device_init(struct simulation *sim,struct device *in)
{
	in->remesh= -1;
	in->newmeshsize= -1;
	in->Jn_y0=NULL;
	in->Jn_y1=NULL;
	in->Jp_y0=NULL;
	in->Jp_y1=NULL;

	in->n_contact_y0=NULL;
	in->n_contact_y1=NULL;
	in->n_contact_x0=NULL;
	in->n_contact_x1=NULL;

	in->Nad= NULL;
	in->G= NULL;
	in->Gn= NULL;
	in->Gp= NULL;
	in->n= NULL;
	in->p= NULL;
	in->dn= NULL;
	in->dndphi= NULL;
	in->dp= NULL;
	in->dpdphi= NULL;
	in->Eg= NULL;
	in->Xi= NULL;
	in->Ev= NULL;
	in->Ec= NULL;
	in->Rfree= NULL;

	in->mun= NULL;
	in->mup= NULL;

	in->Dn= NULL;
	in->Dp= NULL;

	in->epsilonr= NULL;
	in->epsilonr_e0= NULL;

	in->Fn= NULL;
	in->Fp= NULL;
	in->Nc= NULL;
	in->Nv= NULL;

	//Thermal
	in->Tl= NULL;
	in->Te= NULL;
	in->Th= NULL;

	in->Hl= NULL;
	in->H_recombination= NULL;
	in->H_joule= NULL;

	in->He= NULL;
	in->Hh= NULL;

	in->ke= NULL;
	in->kh= NULL;

	in->Fi= NULL;

	in->Fi0_y0 = NULL;
	in->Fi0_y1 = NULL;
	in->Fi0_x0 = NULL;
	in->Fi0_x1 = NULL;

	in->imat= NULL;
	in->imat_epitaxy= NULL;
	in->mask= NULL;

	dim_init_zx_epitaxy(&(in->dim_epitaxy));
	in->mask_epitaxy= NULL;


	in->Jn= NULL;
	in->Jp= NULL;

	in->Jn_diffusion= NULL;
	in->Jn_drift= NULL;

	in->Jp_diffusion= NULL;
	in->Jp_drift= NULL;

	in->Vapplied_y0=NULL;
	in->Vapplied_y1=NULL;
	in->Vapplied_x0=NULL;
	in->Vapplied_x1=NULL;

	in->V_y0= NULL;
	in->V_y1= NULL;
	in->V_x0= NULL;
	in->V_x1= NULL;

	in->t= NULL;
	in->tp= NULL;
	in->kf= NULL;
	in->kd= NULL;
	in->kr= NULL;

	in->Rn= NULL;
	in->Rp= NULL;
	in->Rnet= NULL;

	in->excite_conv= -1;

	in->deltaFln= -1.0;
	in->deltaFlp= -1.0;
	in->deltaFrn= -1.0;
	in->deltaFrp= -1.0;

	in->Rbi_k= NULL;

	in->ex= NULL;
	in->Dex= NULL;
	in->Hex= NULL;

	in->nf_save= NULL;
	in->pf_save= NULL;
	in->nt_save= NULL;
	in->pt_save= NULL;

	in->nfequlib= NULL;
	in->pfequlib= NULL;
	in->ntequlib= NULL;
	in->ptequlib= NULL;

	in->ntb_save= NULL;
	in->ptb_save= NULL;

	in->phi_save= NULL;

	in->xlen= -1.0;
	in->ylen= -1.0;
	in->zlen= -1.0;

	matrix_init(&(in->mx));

//math
	in->max_electrical_itt= -1;
	in->electrical_clamp= -1.0;
	in->max_electrical_itt0= -1;
	in->electrical_clamp0= -1.0;
	in->electrical_error0= -1.0;
	in->math_enable_pos_solver= -1.0;
	in->min_cur_error= -1.0;
	in->pos_max_ittr= -1;
	strcpy(in->solver_name,"");
	strcpy(in->complex_solver_name,"");
	strcpy(in->newton_name,"");

	in->flip_current=1.0;

	in->dt= -1.0;
	in->srh_sim= -1;
	in->go_time= -1;
	in->time= -1.0;
	in->nlast= NULL;
	in->plast= NULL;
	in->ntrapnewton= -1;
	in->ptrapnewton= -1;

	in->stop= -1;
	in->Rshort= -1.0;
	in->onlypos= -1;
	in->odes= -1;
	in->posclamp= -1.0;
	in->wn= NULL;
	in->wp= NULL;

//n traps
	in->nt_all= NULL;
	in->nt= NULL;
	in->ntlast= NULL;
	in->dnt= NULL;
	in->srh_n_r1= NULL;
	in->srh_n_r2= NULL;
	in->srh_n_r3= NULL;
	in->srh_n_r4= NULL;
	in->dsrh_n_r1= NULL;
	in->dsrh_n_r2= NULL;
	in->dsrh_n_r3= NULL;
	in->dsrh_n_r4= NULL;
	in->Fnt= NULL;
	in->tt= NULL;

	in->nt_r1= NULL;
	in->nt_r2= NULL;
	in->nt_r3= NULL;
	in->nt_r4= NULL;
//p traps
	in->pt_all= NULL;
	in->pt= NULL;
	in->ptlast= NULL;
	in->dpt= NULL;
	in->srh_p_r1= NULL;
	in->srh_p_r2= NULL;
	in->srh_p_r3= NULL;
	in->srh_p_r4= NULL;
	in->dsrh_p_r1= NULL;
	in->dsrh_p_r2= NULL;
	in->dsrh_p_r3= NULL;
	in->dsrh_p_r4= NULL;
	in->Fpt= NULL;
	in->tpt= NULL;

	in->pt_r1= NULL;
	in->pt_r2= NULL;
	in->pt_r3= NULL;
	in->pt_r4= NULL;

	in->A= -1.0;
	in->Vol= -1.0;

	in->Rshunt= -1.0;
	in->Rcontact= -1.0;
	in->Rload= -1.0;
	in->L= -1.0;


	in->stop_start= -1;
	in->externalv= -1.0;
	in->Ilast= -1.0;
	in->timedumpcount= -1;
	strcpy(in->simmode,"");
	in->area= -1.0;

	in->nrelax= NULL;
	in->ntrap_to_p= NULL;
	in->prelax= NULL;
	in->ptrap_to_n= NULL;

	in->electrons_y0= NULL;
	in->holes_y0= NULL;
	in->electrons_y1= NULL;
	in->holes_y1= NULL;

	in->electrons_x0= NULL;
	in->holes_x0= NULL;
	in->electrons_x1= NULL;
	in->holes_x1= NULL;

	in->passivate_y0 = NULL;
	in->passivate_y1 = NULL;
	in->passivate_x0 = NULL;
	in->passivate_x1 = NULL;


	in->dumpitdos= -1;


	in->t_big_offset= -1.0;

	in->C= -1.0;
	in->other_layers= -1.0;

	in->kl_in_newton= -1;
	in->config_kl_in_newton= -1;
	in->B= NULL;
	in->xnl_left= -1.0;
	in->xpl_left= -1.0;
	in->stoppoint= -1;
	in->ilast= -1.0;

	in->newton_clever_exit= -1;
	strcpy(in->plot_file,"");

	in->start_stop_time= -1.0;


	in->Is= -1.0;
	in->n_id= -1.0;
	in->Igen= -1.0;

	in->n_orig= NULL;
	in->p_orig= NULL;
	in->n_orig_f= NULL;
	in->p_orig_f= NULL;
	in->n_orig_t= NULL;
	in->p_orig_t= NULL;

	in->Vbi= -1.0;
	in->newton_min_itt= -1;
	in->vbi= -1.0;
	in->avg_gen= -1.0;
	in->dump_energy_slice_xpos= -1;
	in->dump_energy_slice_ypos= -1;
	in->dump_energy_slice_zpos= -1;

	in->pl_intensity= -1.0;
	in->pl_intensity_tot= -1.0;

	in->Rext= -1.0;
	in->Cext= -1.0;
	in->VCext_last= -1.0;
	in->VCext= -1.0;
	in->newton_last_ittr= -1;
	in->phi_mul= -1.0;

	//Newton
	in->newton_dntrap=NULL;
	in->newton_dntrapdntrap=NULL;
	in->newton_dntrapdn=NULL;
	in->newton_dntrapdp=NULL;
	in->newton_dJdtrapn=NULL;
	in->newton_dJpdtrapn=NULL;

	in->newton_dptrapdp=NULL;
	in->newton_dptrapdptrap=NULL;
	in->newton_dptrap=NULL;
	in->newton_dptrapdn=NULL;
	in->newton_dJpdtrapp=NULL;
	in->newton_dJdtrapp=NULL;
	in->newton_dphidntrap=NULL;
	in->newton_dphidptrap=NULL;
	in->newton_ntlast=NULL;
	in->newton_ptlast=NULL;

	in->tm_sun=NULL;
	in->tm_voltage=NULL;
	in->tm_laser=NULL;
	in->tm_time_mesh=NULL;
	in->tm_fs_laser=NULL;
	in->tm_mesh_len=-1;
	in->tm_use_mesh=-1;
	in->tm_mesh_pos=-1;
	in->ncontacts=-1;
	in->active_contact=-1;
	in->newton_only_fill_matrix=FALSE;

	in->dynamic_mesh=-1;

	newton_state_init(&(in->ns));
	newton_state_init(&(in->ns_save));

	//dim_init(&(in->dim_max));
	mesh_obj_init(&(in->mesh_data));
	mesh_obj_init(&(in->mesh_data_save));


	circuit_init(&(in->cir));
	in->circuit_simulation=FALSE;

	in->obj=NULL;
	in->objects=0;
	in->triangles=0;

	shape_init(sim,&(in->big_box));

	//thermal
	heat_init(&(in->thermal));

}


