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
#include "sim.h"
#include "dump.h"
#include "newton_tricks.h"
#include <cal_path.h>
#include <contacts.h>
#include <dump.h>
#include <log.h>
#include <string.h>
#include <stdlib.h>
#include <inp.h>
#include <list.h>
#include <advmath.h>

static int unused __attribute__((unused));

void state_gen_vector(struct simulation *sim,struct device *in)
{
	int i;
	long double *buf;
	int len=0;
	int buf_pos=0;

	char cache_vector[PATH_MAX];
	join_path(3, cache_vector,get_cache_path(sim),in->cache.hash,"vector.dat");

	len+=5;
	len+=17*in->my_epitaxy.electrical_layers;
	len+=2*in->ncontacts;

	buf=malloc(sizeof(long double)*len);

	buf[buf_pos++]=in->zmeshpoints;					//1
	buf[buf_pos++]=in->xmeshpoints;					//2
	buf[buf_pos++]=in->ymeshpoints;					//3
	buf[buf_pos++]=in->srh_bands;					//4
	buf[buf_pos++]=in->my_epitaxy.electrical_layers;	//5

	for (i=0;i<in->my_epitaxy.electrical_layers;i++)
	{
		buf[buf_pos++]=in->dosn[i].config.mu;			//1
		buf[buf_pos++]=in->dosp[i].config.mu;			//2

		buf[buf_pos++]=in->dosn[i].config.doping_start;	//3
		buf[buf_pos++]=in->dosn[i].config.doping_stop;	//4

		buf[buf_pos++]=in->dosn[i].config.srh_sigman;	//5
		buf[buf_pos++]=in->dosn[i].config.srh_sigmap;	//6
		buf[buf_pos++]=in->dosp[i].config.srh_sigman;	//7
		buf[buf_pos++]=in->dosp[i].config.srh_sigmap;	//8

		buf[buf_pos++]=in->dosn[i].config.Nc;			//9
		buf[buf_pos++]=in->dosn[i].config.Nv;			//10

		buf[buf_pos++]=in->dosn[i].config.Eg;			//11
		buf[buf_pos++]=in->dosn[i].config.Xi;			//12
		buf[buf_pos++]=in->dosn[i].config.B;			//13

		buf[buf_pos++]=in->dosn[i].config.Et;			//14
		buf[buf_pos++]=in->dosp[i].config.Et;			//15

		buf[buf_pos++]=in->dosn[i].config.Nt;			//16
		buf[buf_pos++]=in->dosp[i].config.Nt;			//17
	}


	for (i=0;i<in->ncontacts;i++)
	{
		buf[buf_pos++]=in->contacts[i].np;			//1
		buf[buf_pos++]=in->contacts[i].charge_type;	//2
	}


	write_zip_buffer(sim,cache_vector,buf,len);

	free(buf);
}

long double state_load_vector(struct simulation *sim,struct device *in,char *file_name)
{

	int i;
	long double *buf;
	int buf_pos=0;
	long double ret=0.0;

	int len=read_zip_buffer(sim,file_name,&buf);
	if (len==-1)
	{
		return -1;
	}

	
	if (buf[buf_pos++]!=in->zmeshpoints)
	{
		return -1;
	}

	if (buf[buf_pos++]!=in->xmeshpoints)
	{
		return -1;
	}

	if (buf[buf_pos++]!=in->ymeshpoints)
	{
		return -1;
	}

	if (buf[buf_pos++]!=in->srh_bands)
	{
		return -1;
	}

	if (buf[buf_pos++]!=in->my_epitaxy.electrical_layers)
	{
		return -1;
	}

	for (i=0;i<in->my_epitaxy.electrical_layers;i++)
	{
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.mu);	//1
		ret+=log_delta(buf[buf_pos++],in->dosp[i].config.mu);	//2


		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.doping_start);	//3
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.doping_stop);	//4

		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.srh_sigman);	//5
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.srh_sigmap);	//6
		ret+=log_delta(buf[buf_pos++],in->dosp[i].config.srh_sigman);	//7
		ret+=log_delta(buf[buf_pos++],in->dosp[i].config.srh_sigmap);	//8
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.Nc);			//9
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.Nv);			//10
		ret+=fabsl(buf[buf_pos++]-in->dosn[i].config.Eg);				//11
		ret+=fabsl(buf[buf_pos++]-in->dosn[i].config.Xi);				//12
		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.B);			//13
		ret+=fabsl(buf[buf_pos++]-in->dosn[i].config.Et)*1e3;			//14
		ret+=fabsl(buf[buf_pos++]-in->dosp[i].config.Et)*1e3;			//15

		ret+=log_delta(buf[buf_pos++],in->dosn[i].config.Nt);			//16
//		printf("%Le %Le %Le\n",buf[buf_pos-1],in->dosn[i].config.Nt,log_delta(buf[buf_pos-1],in->dosn[i].config.Nt));
		ret+=log_delta(buf[buf_pos++],in->dosp[i].config.Nt);			//17

	}

	for (i=0;i<in->ncontacts;i++)
	{
		ret+=log_delta(buf[buf_pos++],in->contacts[i].np);			//1
		ret+=log_delta(buf[buf_pos++],in->contacts[i].charge_type);	//2
	}


	free(buf);

	return ret;
}

int state_find_vector(struct simulation *sim,struct device *in,char *out)
{

	int i=0;
	char *buffer;
	unsigned int len;
	long l;
	long double min_dv=1e6;
	long double ddevice=0.0;
	long double min_ddevice=1e6;


	struct list files;
	struct inp_file inp;
	char cache_vector[PATH_MAX];
	long double dv;
	char min_file[PATH_MAX];
	char temp_file_name[PATH_MAX];
	int found=FALSE;
	long double error=1e6;
	long double min_error=1e6;

	inp_listdir(sim, get_cache_path(sim),&files);


	for (i=0;i<files.len;i++)
	{
		if ((strcmp(files.names[i],".")!=0)&&(strcmp(files.names[i],"..")!=0))
		{

			join_path(3, cache_vector,get_cache_path(sim),files.names[i],"vector.dat");
			
			ddevice=state_load_vector(sim,in,cache_vector);
			printf("%s %.2Lf\n",files.names[i],ddevice);
			//getchar();
			if (ddevice<1.5)
			{
				if (state_search(sim,in,&dv,files.names[i],temp_file_name,FALSE)==TRUE)
				{
					error=dv+ddevice;
					if (error<min_error)
					{
						min_error=error;
						min_ddevice=ddevice;
						min_dv=dv;
						strcpy(min_file,temp_file_name);
						found=TRUE;
						printf("closest file%s dv=%Lf ddevice=%Lf error=%Lf\n",min_file,dv,ddevice,error);
						//getchar();
					}
				}
			}
		}
	}

	list_free(&files);
	//printf("do load %Lf\n",min_dv);
	//getchar();
	//exit(0);
	if (found==TRUE)
	{

		printf(">>min_file%s dv=%Le ddevice=%Le error=%Le\n",min_file,min_dv,min_ddevice,min_error);
		//getchar();
		if (load_state(sim,in,min_file)==FALSE)
		{
			ewe(sim,"Probem with load\n");
		}
		return TRUE;

	}

	return FALSE;
}
