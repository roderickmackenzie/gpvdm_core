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

/** @file memory.h
@brief allocate 3D memory
*/
#ifndef memory_h
#define memory_h

#include <device.h>

void three_d_copy_gdouble(struct dimensions *dim, gdouble ***dst, gdouble ***src);
void three_d_mul_gdouble(struct dimensions *dim, gdouble ***src, gdouble val);
void malloc_zx_gdouble(struct dimensions *dim, gdouble * (**var));
void free_zx_gdouble(struct dimensions *dim, gdouble **var);
void free_srh_bands(struct dimensions *dim, gdouble **** var);
void malloc_3d_gdouble(struct dimensions *dim, gdouble * (***var));
void free_3d_gdouble(struct dimensions *dim, gdouble ***var);
void malloc_3d_int(struct dimensions *dim, int * (***var));
void free_3d_int(struct dimensions *dim, int ***var);
void malloc_srh_bands(struct dimensions *dim, gdouble * (****var));
void malloc_zx_int(struct dimensions *dim, int * (**var));
void free_zx_int(struct dimensions *dim, int **var);

//3d opps
long double three_d_avg(struct device *in, gdouble ***src);
long double three_d_integrate(struct dimensions *dim, long double ***src);
long double three_d_avg_fabsl(struct device *in, long double ***src);
void three_d_printf(struct dimensions *dim, long double ***src);
void three_d_set_gdouble(struct dimensions *dim, gdouble ***var, gdouble val);
void three_d_sub_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***sub);
void three_d_add_gdouble(struct dimensions *dim, gdouble ***var, gdouble ***add);
void three_d_interpolate_gdouble(long double ***out, long double ***in, struct dimensions *dim_out, struct dimensions *dim_in);
void three_d_quick_dump(char *file_name, long double ***in, struct dimensions *dim);
void three_d_interpolate_srh(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band);
void srh_quick_dump(char *file_name, long double ****in, struct dimensions *dim,int band);
void three_d_interpolate_srh2(long double ****out, long double ****in, struct dimensions *dim_out, struct dimensions *dim_in,int band);
//2d opps
void mem_set_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in);
void mem_add_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in);

void memory_flip_1d_long_double(long double *var,int len);
void memory_flip_1d_int(int *var,int len);

//srh opps
void malloc_srh_bands(struct dimensions *dim, gdouble * (****var));
void srh_copy_gdouble(struct dimensions *dim, gdouble ****dst, gdouble ****src);
#endif
