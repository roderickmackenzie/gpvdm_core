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



/** @file list.h
	@brief Header file for list.c
*/
#ifndef list_h
#define list_h

#include <sim_struct.h>

struct vec_list
{
struct vec *list;
int max;
int length;
double cog_x;
double cog_y;
double max_y;
double min_y;
};

void vec_list_load(struct simulation *sim,struct vec_list* in,char *file_name);
int vec_list_check(struct simulation *sim,struct vec_list* in,struct vec *test);
void vec_list_dump(struct simulation *sim,char *file_name,struct vec_list* in);
void vec_list_init(struct simulation *sim,struct vec_list* in);
void vec_list_add_no_rep(struct simulation *sim,struct vec_list* in,struct vec *test);
void vec_list_add(struct simulation *sim,struct vec_list* in,double one, double two);
int vec_list_get_length(struct simulation *sim,struct vec_list* in);
void vec_list_free(struct simulation *sim,struct vec_list* in);
void vec_list_dump_2d(struct simulation *sim,char *file_name,struct vec_list* in);
void vec_list_cog_cal(struct simulation *sim,struct vec_list* in);
void vec_list_minmax_cal(struct simulation *sim,struct vec_list* in);
void vec_list_remove_bump_up(struct simulation *sim,struct vec_list* in,int start);
void vec_list_remove_bump_down(struct simulation *sim,struct vec_list* in,int start);
#endif

