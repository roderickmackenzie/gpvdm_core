//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
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
//


/** @file i.h
	@brief Header file for i.c
*/
#ifndef i_h
#define i_h
#include "advmath.h"
#include <sim_struct.h>
#include <i_struct.h>

void inter_malloc(struct math_xy* in,int len);
void inter_realloc(struct math_xy* in,int len);
int inter_get_col_n(struct simulation *sim,char *name);
void inter_add_to_hist(struct math_xy* in,gdouble pos,gdouble value);
void inter_init_mesh(struct math_xy* in,int len,gdouble min,gdouble max);
void inter_smooth_range(struct math_xy* out,struct math_xy* in,int points,gdouble x);
gdouble inter_avg_range(struct math_xy* in,gdouble start,gdouble stop);
gdouble inter_array_get_max(gdouble *data,int len);
void inter_div(struct simulation *sim,struct math_xy* one,struct math_xy* two);
void inter_div_long_double(struct math_xy* in,gdouble div);
gdouble inter_get_min_range(struct math_xy* in,gdouble min, gdouble max);
void inter_make_cumulative(struct math_xy* in);
void inter_y_mul_dx(struct math_xy* in);
void inter_add_x(struct math_xy* in,gdouble value);
int inter_sort(struct math_xy* in);
gdouble inter_get_quartile(struct math_xy* in,gdouble value);
void inter_save_seg(struct math_xy* in,char *path,char *name,int seg);
gdouble inter_intergrate(struct math_xy* in);
void inter_to_log_mesh(struct math_xy* out,struct math_xy* in);
void inter_smooth(struct math_xy* out,struct math_xy* in,int points);
gdouble inter_sum_mod(struct math_xy* in);
void inter_set_value(struct math_xy* in,gdouble value);
gdouble inter_get_neg(struct math_xy* in,gdouble x);
gdouble inter_get_noend(struct math_xy* in,gdouble x);
void inter_to_new_mesh(struct math_xy* in,struct math_xy* out);
void inter_swap(struct math_xy* in);
void inter_log_y_m(struct math_xy* in);
gdouble inter_get_min(struct math_xy* in);
gdouble inter_get_fabs_max(struct math_xy* in);
gdouble inter_norm_to_one_range(struct math_xy* in,gdouble min,gdouble max);
void inter_chop(struct math_xy* in,gdouble min, gdouble max);
void inter_save_a(struct math_xy* in,char *path,char *name);
void inter_dump(struct simulation *sim,struct math_xy* in);
void inter_purge_zero(struct math_xy* in);
void inter_append(struct math_xy* in,gdouble x,gdouble y);
void inter_init(struct simulation *sim,struct math_xy* in);
void inter_sub_long_double(struct math_xy* in,gdouble value);
void inter_sub(struct simulation *sim,struct math_xy* one,struct math_xy* two);
gdouble inter_sum(struct math_xy* in);
void inter_copy(struct math_xy* in,struct math_xy* orig,int alloc);
int inter_get_col(char *file);
void inter_load_by_col(struct simulation *sim,struct math_xy* in,char *name,int col);
gdouble inter_get_diff(char *out_path,struct math_xy* one,struct math_xy* two,gdouble start,gdouble stop,struct math_xy* mull);
void inter_pow(struct math_xy* in,gdouble p);
gdouble inter_get_raw(gdouble *x,gdouble *data,int len,gdouble pos);
gdouble inter_norm(struct math_xy* in,gdouble mul);
void inter_log_y(struct math_xy* in);
void inter_mul(struct math_xy* in,gdouble mul);
void inter_log_x(struct math_xy* in);
void inter_save(struct math_xy* in,char *name);
int inter_load(struct simulation *sim,struct math_xy* in,char *name);
gdouble inter_get_hard(struct math_xy* in,gdouble x);
gdouble inter_get(struct math_xy* in,gdouble x);
void inter_print(struct math_xy* in);
void inter_free(struct math_xy* in);
void inter_rescale(struct math_xy* in,gdouble xmul, gdouble ymul);
void inter_mod(struct math_xy* in);
void inter_add(struct math_xy* out,struct math_xy* in);
void inter_norm_area(struct math_xy* in,gdouble mul);
gdouble inter_get_max(struct math_xy* in);
gdouble inter_get_max_range(struct math_xy* in,int start, int stop);
void inter_add_long_double(struct math_xy* in,gdouble value);
gdouble inter_intergrate_lim(struct math_xy* in,gdouble from, gdouble to);
void inter_deriv(struct math_xy* out,struct math_xy* in);
void inter_import_array(struct math_xy* in,gdouble *x,gdouble *y,int len,int alloc);
gdouble inter_avg(struct math_xy* in);
void inter_convolve(struct math_xy* one,struct math_xy* two);
void inter_save_backup(struct math_xy* in,char *name,int backup);
void inter_dft(gdouble *real,gdouble *imag,struct math_xy* in,gdouble fx);
int inter_get_max_pos(struct math_xy* in);
int inter_search_pos(struct math_xy* in,gdouble x);
void inter_join_bins(struct math_xy* in,gdouble delta);
void inter_reset(struct math_xy* in);
void inter_find_peaks(struct math_xy* out,struct math_xy* in,int find_max);
void inter_sin(struct math_xy *in,gdouble mag,gdouble fx,gdouble delta);
void inter_purge_x_zero(struct math_xy* in);
int inter_search_token(struct simulation *sim,long double *value,char *token,char *name);
int inter_get_min_pos(struct math_xy* in);
void math_xy_get_left_right_start(struct math_xy* in,int *left,int *right, long double fraction);
#endif
