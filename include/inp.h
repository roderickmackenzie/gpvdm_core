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

/** @file inp.h
@brief Code to read input files.
*/

#ifndef inp_h
#define inp_h
#include "advmath.h"
#include "inp_struct.h"
#include <sim_struct.h>
#include "list_struct.h"

void inp_replace_double(struct simulation *sim,struct inp_file *in,char *token, double value);
int inp_get_array_len(struct simulation *sim,struct inp_file *in,char *token);
int inp_get_array(struct simulation *sim,char ** out,struct inp_file *in,char *token);
int inp_save(struct simulation *sim,struct inp_file *in);
void inp_init(struct simulation *sim,struct inp_file *in);
int inp_aes_load(struct inp_file *in,char *path,char *file,char *key);
int inp_load(struct simulation *sim,struct inp_file *in,char *file);
void inp_free(struct simulation *sim,struct inp_file *in);
int inp_search_gdouble(struct simulation *sim,struct inp_file *in,gdouble* out,char* token);
int inp_search_double(struct simulation *sim,struct inp_file *in,double* out,char* token);
int inp_search_float(struct simulation *sim,struct inp_file *in,float* out,char* token);
void inp_search_int(struct simulation *sim,struct inp_file *in,int* out,char* token);
void inp_search_string(struct simulation *sim,struct inp_file *in,char* out,char* token);
int inp_search(struct simulation *sim,char* out,struct inp_file *in,char *token);
void inp_check(struct simulation *sim,struct inp_file *in,double ver);
int inp_read_buffer(struct simulation *sim,char **buf, long *len,char *full_file_name);
void inp_reset_read(struct simulation *sim,struct inp_file *in);
char* inp_get_string(struct simulation *sim,struct inp_file *in);
char* inp_search_part(struct simulation *sim,struct inp_file *in,char *token);
int inp_load_from_path(struct simulation *sim,struct inp_file *in,char *path,char *file);
void inp_replace(struct simulation *sim,struct inp_file *in,char *token, char *text);
int inp_search_pos(struct simulation *sim,struct inp_file *in,char *token);
int inp_search_english(struct simulation *sim,struct inp_file *in,char *token);
int inp_isfile(struct simulation *sim,char *full_file_name);
int zip_is_in_archive(char *full_file_name);
int zip_write_buffer(struct simulation *sim,char *full_file_name,char *buffer, int len);

void inp_listdir(struct simulation *sim, char *dir_name,struct list *out);
int guess_whole_sim_name(struct simulation *sim,char *ret,char *dir_name,char* search_name);
int search_for_token(struct simulation *sim,char *ret,char *dir_name,char* token,char *search_value);

void inp_replace_offset(struct simulation *sim,struct inp_file *in,char *token, char *text,int offset);
void inp_replace_double_offset(struct simulation *sim,struct inp_file *in,char *token, double value,int offset);
void inp_search_double_offset(struct simulation *sim,struct inp_file *in,double* out,char* token,int offset);
int inp_get_array_gdouble(struct simulation *sim,long double * out,struct inp_file *in,char *token);
int inp_count_hash_tags(struct simulation *sim,struct inp_file *in);
void inp_replace_double_offset_col(struct simulation *sim,struct inp_file *in,char *token, double value,int offset,int col);
int inp_test_end_of_data(char *line);
void inp_load_from_buffer(struct simulation *sim,struct inp_file *in,char *file,char *buffer,int len);
void inp_search_long_long(struct simulation *sim,struct inp_file *in,long long* out,char* token);

int get_line(char *out,char *data,int len,int *pos);
void inp_free_no_save(struct simulation *sim,struct inp_file *in);
#endif
