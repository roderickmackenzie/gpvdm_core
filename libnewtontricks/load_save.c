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

/** @file ntricks.c
	@brief A collection of helper functions for the newton solver.
*/


#include <stdio.h>
#include <exp.h>
#include "sim.h"
#include "dump.h"
#include "newton_tricks.h"
#include "gui_hooks.h"
#include <plot.h>
#include <cal_path.h>
#include <thermal.h>
#include <contacts.h>
#include <dump.h>
#include <log.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <inp.h>

static int unused __attribute__((unused));

void state_cache_enable(struct simulation *sim,struct device *in)
{
	char temp[200];
	struct inp_file inp;
	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_gpvdm_local_path(sim),"cache.inp");

	inp_search_string(sim,&inp,temp,"#cache_enabled");
	in->cache.enabled=english_to_bin(sim,temp);

	inp_free(sim,&inp);

	if (strcmp(in->newton_name,"newton_simple")==0)
	{
		in->cache.enabled=FALSE;
	}

}
void state_cache_init(struct simulation *sim,struct device *in)
{
	char cache_path[PATH_MAX];

	hash_dir(sim,in->cache.hash);

	gpvdm_mkdir(get_cache_path(sim));

	join_path(2, cache_path,get_cache_path(sim),in->cache.hash);
	gpvdm_mkdir(cache_path);

	state_gen_vector(sim,in);

	in->cache.enabled=FALSE;
}

void state_new_name(struct simulation *sim,struct device *in,char *out_file)
{
	int i;
	FILE *state;

	char fname[100];
	char path[PATH_MAX];

	i=0;
	while(1)
	{
		sprintf(fname,"state%d.dat",i);
		join_path(3, path,get_cache_path(sim),in->cache.hash,fname);
		state=fopen(path,"r");
		if (state!=NULL)
		{
			fclose(state);
		}else
		{
			strcpy(out_file,fname);
			return;
		}

		i++;
	}
}


void save_state(struct simulation *sim,struct device *in)
{
	if (in->cache.enabled==FALSE)
	{
		return;
	}

	if (in->last_error>1e-5)
	{
		return;
	}

	int pos=0;
	int len=0;
	int found=FALSE;
	char cache_file[PATH_MAX];
	char cache_index[PATH_MAX];

	long double *buf;
	FILE *file;

	char fname[100];
	if (state_search(sim,in,NULL,in->cache.hash,cache_file,TRUE)==TRUE)
	{
		printf("FOUND!\n");
		return;
	}

	int i=0;
	int z=0;
	int x=0;
	int y=0;

	int band;
	FILE *state;

	state_new_name(sim,in,fname);

	join_path(3, cache_file,get_cache_path(sim),in->cache.hash,fname);
	join_path(3, cache_index,get_cache_path(sim),in->cache.hash,"index.dat");

	pos=0;
	len=in->ncontacts;

	buf=malloc(sizeof(long double)*len);
	for (i=0;i<in->ncontacts;i++)
	{
		buf[pos++]=in->contacts[i].voltage;
	}

	file=fopen(cache_index,"ab");
	fwrite ( (char*)buf, len*sizeof(long double),1,file);
	fclose(file);
	free(buf);

	pos=0;
	len=in->xmeshpoints*in->ymeshpoints*in->zmeshpoints*3;
	len+=in->xmeshpoints*in->ymeshpoints*in->zmeshpoints*in->srh_bands*2;
	len+=in->ncontacts;
	len+=1;	//in->last_ittr
	len+=1; //in->last_error

	buf=malloc(sizeof(long double)*len);

	for (i=0;i<in->ncontacts;i++)
	{
		buf[pos++]=in->contacts[i].voltage;
	}

	buf[pos++]=(long double)in->last_ittr;
	buf[pos++]=(long double)in->last_error;

	for (z=0;z<in->zmeshpoints;z++)
	{
		for (x=0;x<in->xmeshpoints;x++)
		{
			for (y=0;y<in->ymeshpoints;y++)
			{
				buf[pos++]=in->phi[z][x][y];
				buf[pos++]=in->x[z][x][y];
				buf[pos++]=in->xp[z][x][y];


				for (band=0;band<in->srh_bands;band++)
				{
					buf[pos++]=in->xt[z][x][y][band];
					buf[pos++]=in->xpt[z][x][y][band];
				}

			}
		}
	}

	//printf("%s\n",cache_file);
	write_zip_buffer(sim,cache_file,buf,len);
//	file=fopen(cache_file,"wb");
//	fwrite ( (char*)buf, len*sizeof(long double),1,file);
//	fclose(file);
	free(buf);
}


int state_search_and_load(struct simulation *sim,struct device *in)
{
	if (in->cache.enabled==FALSE)
	{
		return FALSE;
	}

	int found=FALSE;
	char cache_file[PATH_MAX];
	if (state_search(sim,in,NULL,in->cache.hash,cache_file,TRUE)==TRUE)
	{
		if (load_state(sim,in,cache_file)==FALSE)
		{
			ewe(sim,"Probem with load\n");
		}
		return TRUE;
	}

	/*
	if (non_exact_answer==TRUE)
	{
		if (load_state(sim,in,cache_file)==FALSE)
		{
			ewe(sim,"Probem with load\n");
		}
		return FALSE;
	}*/

return FALSE;
}

struct state
{
	long double Vc;		//Contact voltage
	long double np;		//Contact carrier density
	int charge_type;			//Contact type
};

int state_search(struct simulation *sim,struct device *in,long double *ret_error,char *hash_dir,char *file_name,int actual)
{
	if (in->cache.enabled==FALSE)
	{
		return FALSE;
	}

	int x;
	int y;
	int z;
	int i;
	int band;
	long double tot_error;
	long double V;
	FILE *state;
	char path[PATH_MAX];
	char full_path[PATH_MAX];
	int test=FALSE;
	DIR *d;
	struct dirent *dir;
	char cache_index[PATH_MAX];
	long double min_error=1e6;
	int min_pos=0;
	struct state states[100];
	char cache_file[100];
	long double *buf;
	FILE *file;
    struct stat st;
	int item;
	int items;
	int pos=0;
	int ret=FALSE;

	//printf("\n");

	//getchar();
	join_path(3, cache_index,get_cache_path(sim),hash_dir,"index.dat");

    if (stat(cache_index,&st)!=0)
	{
		return FALSE;
	}

	int len=st.st_size;
	items=len/(sizeof(long double)*in->ncontacts);
	
	buf=malloc(len);

	file=fopen(cache_index,"rb");
	fread((long double*)buf, len, 1, file);
	fclose(file);


	for (item=0;item<items;item++)
	{
	//printf("items=%d %d\n",item,items);

		pos=item*in->ncontacts;
		tot_error=0.0;

		for (i=0;i<in->ncontacts;i++)
		{	
			V=buf[pos+i];
			if (actual==TRUE)
			{
				tot_error+=fabsl(in->contacts[i].voltage-V);
			}else
			{
				tot_error+=fabsl(in->contacts[i].voltage_want-V);
			}
			//printf("%.20Le %.20Le %Le\n",in->contacts[i].voltage,V,in->contacts[i].voltage-V);
			//getchar();
		}


//		printf("%d %Le\n",item,tot_error);

		if (tot_error==0.0)
		{
			sprintf(cache_file,"state%d.dat",item);
			join_path(3, full_path,get_cache_path(sim),hash_dir,cache_file);
			strcpy(file_name,full_path);
			//printf("%s\n",file_name);
			//getchar();
			return TRUE;
		}

		if (ret_error!=NULL)
		{
			//printf("%Le\n",tot_error);

			if (tot_error<min_error)
			{
				min_error=tot_error;
				*ret_error=min_error;
				sprintf(cache_file,"state%d.dat",item);
				join_path(3, full_path,get_cache_path(sim),hash_dir,cache_file);
				strcpy(file_name,full_path);
				ret=TRUE;
			}
		}

	}

	free(buf);

//getchar();
return ret;
}

int load_state(struct simulation *sim,struct device *in,char *file_name)
{

	//printf_log(sim,"Load state %s\n",file_name);
	int i=0;
	int z=0;
	int x=0;
	int y=0;

	int band=0;
	long double V;
	gdouble vtest;
	int len=0;

	long double *buf;
	len=read_zip_buffer(sim,file_name,&buf);

	if (len==-1)
	{
		getchar();
		ewe(sim,"load_state file not found %s\n",file_name);		
	}

	int pos=0;

	for (i=0;i<in->ncontacts;i++)
	{
		in->contacts[i].voltage=buf[pos++];
	}

	in->last_ittr=(int)buf[pos++];
	in->last_error=(long double)buf[pos++];

	contacts_update(sim,in);

	//printf_log(sim,"Loading state\n");

	for (z=0;z<in->zmeshpoints;z++)
	{
		for (x=0;x<in->xmeshpoints;x++)
		{
			for (y=0;y<in->ymeshpoints;y++)
			{
				in->phi[z][x][y]=buf[pos++];
				in->x[z][x][y]=buf[pos++];
				in->xp[z][x][y]=buf[pos++];


				for (band=0;band<in->srh_bands;band++)
				{
					in->xt[z][x][y][band]=buf[pos++];
					in->xpt[z][x][y][band]=buf[pos++];
				}

			}

			//printf("here\n");
			update_y_array(sim,in,z,x);
		}
	}

	free(buf);

return TRUE;

}


