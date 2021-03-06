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

/** @file dump_energy_slice.c
@brief dump an energy slice, code needs updating
*/

#include <string.h>
#include <stdio.h>
#include <dos.h>
#include "sim.h"
#include "dump.h"
#include "dat_file.h"

void dump_energy_slice(struct simulation *sim,char *out_dir,struct device *in)
{
struct dat_file buf;
buffer_init(&buf);
char name[200];
char outpath[200];
int band=0;
int x=in->dump_energy_slice_xpos;
int y=in->dump_energy_slice_ypos;
int z=in->dump_energy_slice_zpos;
struct dimensions *dim=&in->ns.dim;

if ((y>=dim->ylen)||(x>=dim->xlen)||(z>=dim->zlen))
{
	return;
}

struct math_xy dump1;
inter_init(sim,&dump1);

struct math_xy dump2;
inter_init(sim,&dump2);

struct math_xy dump3;
inter_init(sim,&dump3);

struct math_xy dump4;
inter_init(sim,&dump4);

struct math_xy dump5;
inter_init(sim,&dump5);

struct math_xy dump6;
inter_init(sim,&dump6);

struct math_xy dump9;
inter_init(sim,&dump9);

struct math_xy dump10;
inter_init(sim,&dump10);

struct math_xy dump_nt;
inter_init(sim,&dump_nt);

struct math_xy dump_pt;
inter_init(sim,&dump_pt);

int mat=in->imat[z][x][y];
for (band=0;band<dim->srh_bands;band++)
{
	inter_append(&dump1,get_dos_E_n(in,band,mat),in->n[z][x][y]*in->srh_n_r1[z][x][y][band]-in->srh_n_r2[z][x][y][band]);
	inter_append(&dump2,get_dos_E_p(in,band,mat),in->p[z][x][y]*in->srh_p_r1[z][x][y][band]-in->srh_p_r2[z][x][y][band]);
	inter_append(&dump3,get_dos_E_n(in,band,mat),in->nt[z][x][y][band]-in->ntb_save[z][x][y][band]);
	inter_append(&dump4,get_dos_E_p(in,band,mat),in->pt[z][x][y][band]-in->ptb_save[z][x][y][band]);
	inter_append(&dump5,get_dos_E_n(in,band,mat),in->p[z][x][y]*in->srh_n_r3[z][x][y][band]-in->srh_n_r4[z][x][y][band]);
	inter_append(&dump6,get_dos_E_p(in,band,mat),in->n[z][x][y]*in->srh_p_r3[z][x][y][band]-in->srh_p_r4[z][x][y][band]);
	inter_append(&dump9,get_dos_E_n(in,band,mat),in->Fnt[z][x][y][band]);
	inter_append(&dump10,get_dos_E_p(in,band,mat),in->Fpt[z][x][y][band]);
	inter_append(&dump_nt,get_dos_E_n(in,band,mat),in->nt[z][x][y][band]);
	inter_append(&dump_pt,get_dos_E_p(in,band,mat),in->pt[z][x][y][band]);

}

sprintf(outpath,"%senergy_slice_nt_cap_%s",out_dir,".dat");
//inter_save(&dump1,outpath);
inter_free(&dump1);

sprintf(outpath,"%senergy_slice_pt_cap_%s",out_dir,".dat");
//inter_save(&dump2,outpath);
inter_free(&dump2);


sprintf(outpath,"%senergy_slice_nt_delta_%s",out_dir,".dat");
//inter_save(&dump3,outpath);
inter_free(&dump3);

sprintf(outpath,"%senergy_slice_pt_delta_%s",out_dir,".dat");
//inter_save(&dump4,outpath);
inter_free(&dump4);


sprintf(outpath,"%senergy_slice_nt_recom_%s",out_dir,".dat");
//inter_save(&dump5,outpath);
inter_free(&dump5);

sprintf(outpath,"%senergy_slice_pt_recom_%s",out_dir,".dat");
//inter_save(&dump6,outpath);
inter_free(&dump6);

sprintf(outpath,"%senergy_slice_fn_%s",out_dir,".dat");
//inter_save(&dump9,outpath);
inter_free(&dump9);

sprintf(outpath,"%senergy_slice_fp_%s",out_dir,".dat");
//inter_save(&dump10,outpath);
inter_free(&dump10);

buffer_malloc(&buf);
sprintf(name,"energy_slice_nt.dat");
buf.y_mul=1.0;
buf.x_mul=1.0;
strcpy(buf.title,"Energy - trap ocupation");
strcpy(buf.type,"xy");
strcpy(buf.x_label,"Energy");
strcpy(buf.data_label,"Ocupation");
strcpy(buf.x_units,"eV");
strcpy(buf.data_units,"m^{-3} eV^{-1}");
buf.logscale_x=0;
buf.logscale_y=1;
buf.x=1;
buf.y=dump_nt.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,dump_nt.x, dump_nt.data, dump_nt.len);
buffer_dump_path(sim,out_dir,name,&buf);
buffer_free(&buf);
inter_free(&dump_nt);

buffer_malloc(&buf);
sprintf(name,"energy_slice_pt.dat");
buf.y_mul=1.0;
buf.x_mul=1.0;
strcpy(buf.title,"Energy - trap ocupation");
strcpy(buf.type,"xy");
strcpy(buf.x_label,"Energy");
strcpy(buf.data_label,"Ocupation");
strcpy(buf.x_units,"eV");
strcpy(buf.data_units,"m^{-3} eV^{-1}");
buf.logscale_x=0;
buf.logscale_y=1;
buf.x=1;
buf.y=dump_pt.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,dump_pt.x, dump_pt.data, dump_pt.len);
buffer_dump_path(sim,out_dir,name,&buf);
buffer_free(&buf);
inter_free(&dump_pt);

}
