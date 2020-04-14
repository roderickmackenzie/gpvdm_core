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

/** @file triangle.h
@brief triangle structure
*/
#ifndef triangle_h
#define triangle_h

#include <vec.h>

struct triangle
{
	struct vec xy0;
	struct vec xy1;
	struct vec xy2;
	int object_uid;

	int object_type;

	//pre calculation
	struct vec edge1;
	struct vec edge2;

	int obj_left;
	int obj_right;

	//not used
	int flag;
};

struct triangles
{
	struct triangle *data;
	int edges_calculated;
	int max_len;
	int len;
};

void triangle_init(struct triangle *tri);
void triangle_print(struct triangle *in);
void triangle_cog(struct vec *out,struct triangle *in);
void triangle_norm(struct vec *ret,struct triangle *my_obj);
void triangle_dump(char *file_name,struct triangle *tri);
int triangle_vec_within (struct triangle *tri,struct vec *pt);
double triangle_get_y_from_xz(struct triangle *tri,double x, double z);

double triangle_get_min_y(struct triangle* tri);
double triangle_get_max_y(struct triangle* tri);
void ray_tri_flip_y_axis(struct triangle* tri,double y);
void ray_tri_sub_y(struct triangle* tri,double y);
void ray_tri_add_y(struct triangle* tri,double y);
void triangle_cal_edges(struct triangle* tri);
double triangle_cal_area(struct triangle* tri);

#endif
