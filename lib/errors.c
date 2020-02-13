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

/** @file exit.c
	@brief Exit while sending sane data to the gui, also handle crashes while fitting.
*/


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "log.h"
#include "gpvdm_const.h"
//#include "dump_ctrl.h"
//#include "dos_types.h"
#include "gui_hooks.h"
#include "lang.h"
#include <signal.h>
#include <server.h>
#include <unistd.h>

static char lock_name[100];
static char lock_data[100];

void errors_add(struct simulation *sim, const char *format, ...)
{

	char temp[1000];
	char temp2[1000];
	int len=0;
	int next_size=0;
	va_list args;
	va_start(args, format);
	vsprintf(temp,format, args);

	sprintf(temp2,"error:%s\n",temp);
	len=strlen(temp2);

	ewe(sim,temp2);

	next_size=sim->error_log_size+len;
	if (next_size>sim->error_log_size_max)
	{
		sim->error_log_size_max+=1024;
		sim->error_log=realloc(sim->error_log,sizeof(char)*sim->error_log_size_max);
	}

	sim->error_log_size=next_size;
	strcat(sim->error_log,temp2);
	sim->errors_logged++;

}

void errors_dump(struct simulation *sim)
{
	printf("%s\n",sim->error_log);
}

int is_errors(struct simulation *sim)
{
	if (sim->errors_logged==0)
	{
		return -1;
	}

	return 0;
}

void errors_init(struct simulation *sim)
{
	sim->error_log_size=0;
	sim->error_log_size_max=1024;
	sim->error_log=malloc(sizeof(char)*sim->error_log_size_max);
	strcpy(sim->error_log,"");
	sim->errors_logged=0;

}

void errors_free(struct simulation *sim)
{
	sim->error_log_size=0;
	sim->error_log_size_max=0;
	free(sim->error_log);
	sim->errors_logged=0;
}

void set_ewe_lock_file(char *lockname,char *data)
{
strcpy(lock_name,lockname);
strcpy(lock_data,data);
}

int ewe( struct simulation *sim, const char *format, ...)
{
	FILE* out;
	char temp[1000];
	char temp2[1000];
	va_list args;
	va_start(args, format);
	vsprintf(temp,format, args);

	sprintf(temp2,"error:%s",temp);

	printf_log(sim, "%s\n",temp2);

	va_end(args);





	gui_send_finished_to_gui(sim);


	if (strcmp(lock_name,"")!=0)
	{
		out=fopen(lock_name,"w");
		fprintf(out,"%s",lock_data);
		fclose(out);
	}
	//printf_log(sim,"Raising segmentation fault\n");
	//raise(SIGSEGV);
exit(1);
}


