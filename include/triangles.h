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

/** @file triangles.h
@brief Functions to manipulate lots of triangles
*/
#ifndef triangles_io_h
#define triangles_io_h

#include <vec.h>
#include <triangle.h>
void triangle_load_from_file(struct simulation *sim,struct triangles *in,char *file_name);
void triangle_print(struct triangle *in);
void triangles_print(struct triangles *in);
void triangles_free(struct triangles *in);
void triangles_cpy(struct triangles *out,struct triangles *in);
void triangles_find_min(struct vec *out,struct triangles *in);
void triangles_find_max(struct vec *out,struct triangles *in);
void triangles_sub_vec(struct triangles *in,struct vec *v);
void triangles_add_vec(struct triangles *in,struct vec *v);
void triangles_div_vec(struct triangles *in,struct vec *v);
void triangles_mul_vec(struct triangles *in,struct vec *v);
void triangles_cal_edges(struct triangles *in);
void triangles_init(struct triangles *tri);
void triangles_malloc(struct triangles *tri);
void triangles_save(char *file_name,struct triangles *in);
void triangles_add_triangle(struct triangles *obj, double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,int uid,int object_type);
void triangles_set_object_type(struct triangles *in,int object_type);
double triangles_interpolate(struct triangles *in,struct vec *p);

//flags
void triangles_set_flag(struct triangles *obj,int flag);
void triangles_unset_flags_for_flat(struct triangles *obj);
double triangle_Ra(struct simulation *sim,struct triangles *obj);

//Roughness
double triangle_Rq(struct simulation *sim,struct triangles *obj);
#endif
