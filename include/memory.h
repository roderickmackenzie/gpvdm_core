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

/** @file memory.h
@brief allocate 3D memory
*/
#ifndef memory_h
#define memory_h

#include <device.h>

void three_d_copy_gdouble(struct dimensions *dim, gdouble ***dst, gdouble ***src);
void zxy_mul_gdouble(struct dimensions *dim, gdouble ***src, gdouble val);
void zxy_div_gdouble(struct dimensions *dim, gdouble ***src, gdouble val);
void zxy_long_double_mul_by_zxy_long_double(struct dimensions *dim, gdouble ***a, gdouble ***b);
void zxy_long_double_div_by_zxy_long_double(struct dimensions *dim, gdouble ***a, gdouble ***b);

void malloc_zx_gdouble(struct dimensions *dim, gdouble * (**var));
void free_zx_gdouble(struct dimensions *dim, gdouble * (**var));

void free_srh_bands(struct dimensions *dim, gdouble * (**** in_var));
void malloc_srh_bands(struct dimensions *dim, gdouble * (****var));


//zxy_int
void malloc_3d_int(struct dimensions *dim, int * (***var));
void free_3d_int(struct dimensions *dim, int ***var);
void dump_zxy_int(struct dimensions *dim, int ***var);

void malloc_zx_int(struct dimensions *dim, int * (**var));
void free_zx_int(struct dimensions *dim, int *(**in_var));

//3d opps
long double three_d_avg(struct device *in, gdouble ***src);
long double three_d_avg_raw(struct device *in, long double ***src);
long double three_d_integrate(struct dimensions *dim, long double ***src);
long double three_d_avg_fabsl(struct device *in, long double ***src);
void three_d_printf(struct dimensions *dim, long double ***src);
void three_d_sub_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***sub);
void three_d_add_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***add);
void three_d_interpolate_gdouble(long double ***out, long double ***in, struct dimensions *dim_out, struct dimensions *dim_in);
void three_d_quick_dump(char *file_name, long double ***in, struct dimensions *dim);
void three_d_interpolate_srh(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band);
void srh_quick_dump(char *file_name, long double ****in, struct dimensions *dim,int band);
void three_d_interpolate_srh2(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band);
long double zxy_min_gdouble(struct dimensions *dim, gdouble ***var);
long double zxy_max_gdouble(struct dimensions *dim, gdouble ***var);
long double zxy_sum_gdouble(struct dimensions *dim, long double ***src);
long double zx_y_max_gdouble(struct dimensions *dim, gdouble ***var,int y);

//zy_long_double
void malloc_zy_long_double(struct dimensions *dim, long double * (**var));
void free_zy_long_double(struct dimensions *dim, long double * (**in_var));
void malloc_zy_int(struct dimensions *dim, int * (**var));
void free_zy_int(struct dimensions *dim, int *(**in_var));

//zxy_long_double
void malloc_zxy_gdouble(struct dimensions *dim, gdouble * (***var));
void free_zxy_gdouble(struct dimensions *dim, gdouble * (***in_var));
void zxy_load_long_double(struct simulation *sim, struct dimensions *dim,long double * *** data,char *file_name);
void zx_y_quick_dump(char *file_name, long double ***in, struct dimensions *dim);
void zxy_set_gdouble(struct dimensions *dim, gdouble ***var, gdouble val);
void flip_zxy_long_double_y(struct simulation *sim, struct dimensions *dim,long double *** data);
long double interpolate_zxy_long_double(struct dimensions *dim, long double ***data,int z, int x, long double y_in);

//zx_epitaxy_int
void malloc_zx_epitaxy_int(struct dim_zx_epitaxy *dim, int * (***var));
void free_zx_epitaxy_int(struct dim_zx_epitaxy *dim, int ***var);
void dump_zx_epitaxy_int(struct dim_zx_epitaxy *dim, int ***var);

//zxy_long_double_complex
void malloc_zxy_long_double_complex(struct dimensions *dim, long double complex * (***var));
void free_zxy_long_double_complex(struct dimensions *dim, long double complex * (***in_var));

//light_zxyl_long_double
void malloc_light_zxyl_long_double(struct dim_light *dim, long double * (****var));
void free_light_zxyl_long_double(struct dim_light *dim, long double * (****in_var));
void flip_light_zxyl_long_double_y(struct simulation *sim, struct dim_light *dim,long double **** data);
void div_light_zxyl_long_double(struct dim_light *dim, long double ****data,long double val);
void memset_light_zxyl_long_double(struct dim_light *dim, long double ****data,int val);
void memset_light_zxyl_long_double_y(struct dim_light *dim, long double ****data,int z, int x, int l,long double val);

//light_zxy_long_double
void malloc_light_zxy_long_double(struct dim_light *dim, long double * (***var));
void free_light_zxy_long_double(struct dim_light *dim, long double * (***in_var));
void flip_light_zxy_long_double_y(struct simulation *sim, struct dim_light *dim,long double *** data);
void memset_light_zxy_long_double(struct dim_light *dim, long double ***data,int val);
void div_light_zxy_long_double(struct dim_light *dim, long double ***data,long double val);
long double interpolate_light_zxy_long_double(struct dim_light *dim, long double ***data,int z, int x, long double y_in);
long double interpolate_light_zxy_long_double_intergral(struct dim_light *dim, long double ***data,int z, int x, long double y_start,long double y_stop);

//light_l_long_double
void malloc_light_l_long_double(struct dim_light *dim, long double * (*var));
void free_light_l_long_double(struct dim_light *dim, long double * (*in_var));

//light_zxyl_long_double_complex
void malloc_light_zxyl_long_double_complex(struct dim_light *dim, long double complex * (****var));
void free_light_zxyl_long_double_complex(struct dim_light *dim, long double complex * (****in_var));

//light_zxy_p_object
void malloc_light_zxy_p_object(struct dim_light *dim, struct object * (****var));
void free_light_zxy_p_object(struct dim_light *dim, struct object * (****in_var));

// heat_zxy_long_double
void malloc_heat_zxy_long_double(struct dim_heat *dim, long double * (***var));
void free_heat_zxy_long_double(struct dim_heat *dim, long double * (***in_var));
void flip_heat_zxy_long_double_y(struct simulation *sim, struct dim_heat *dim,long double *** data);
void memset_heat_zxy_long_double(struct dim_heat *dim, long double ***data,int val);
void div_heat_zxy_long_double(struct dim_heat *dim, long double ***data,long double val);
long double interpolate_heat_zxy_long_double(struct dim_heat *dim, long double ***data,int z, int x, long double y_in);
long double intergrate_heat_zxy_long_double(struct dim_heat *dim, long double ***data);
long double avg_heat_zxy_long_double(struct dim_heat *dim, long double ***data);

//heat_zxy_p_object
void malloc_heat_zxy_p_object(struct dim_heat *dim, struct object * (****var));
void free_heat_zxy_p_object(struct dim_heat *dim, struct object * (****in_var));

//2d opps
void mem_set_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in);
void mem_add_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in);
void zx_copy_gdouble(struct dimensions *dim, gdouble **dst, gdouble **src);

void memory_flip_1d_long_double(long double *var,int len);
void memory_flip_1d_int(int *var,int len);

//srh opps
void malloc_srh_bands(struct dimensions *dim, gdouble * (****var));
void srh_copy_gdouble(struct dimensions *dim, gdouble ****dst, gdouble ****src);

//matrix
void matrix_init(struct matrix *mx);
void matrix_dump(struct simulation *sim,struct matrix *mx);
void matrix_malloc(struct simulation *sim,struct matrix *mx);
void matrix_free(struct simulation *sim,struct matrix *mx);
void matrix_realloc(struct simulation *sim,struct matrix *mx);
int matrix_solve(struct simulation *sim,struct matrix *mx);
void matrix_cache_reset(struct simulation *sim,struct matrix *mx);
void matrix_save(struct simulation *sim,struct matrix *mx);
int matrix_load(struct simulation *sim,struct matrix *mx);
long double matrix_cal_error(struct simulation *sim,struct matrix *mx);
void matrix_zero_b(struct simulation *sim,struct matrix *mx);
void matrix_dump_b(struct simulation *sim,struct matrix *mx);
void matrix_dump_J(struct simulation *sim,struct matrix *mx);
void matrix_add_nz_item(struct simulation *sim,struct matrix *mx,int x,int y,long double val);
void matrix_convert_J_to_sparse(struct simulation *sim,struct matrix *mx);
void matrix_stats(struct simulation *sim,struct matrix *mx);
//raw memory opps
int search(long double *x,int N,long double find);

#endif
