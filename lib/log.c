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

/** @file log.c
	@brief Deal with output to the log file, can either be disk/screen or both.  UTF8 supported/HTML supported and can stream to GUI under Windows/Linux.
*/



#include <stdio.h>
#include <stdarg.h>
#include "log.h"
#include <colors.h>
#include <time.h>
#include <stdlib.h>
#include <util.h>
#include <cal_path.h>
#include <string.h>
#include <gpvdm_const.h>
#include <lang.h>
#include <dump.h>
#include <wchar.h>
#include <color.h>


void log_time_stamp(struct simulation *sim)
{
	time_t t=0;
	t=time(NULL);

	struct tm tm = *localtime(&t);

	printf_log(sim,"time: %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
}

void log_clear(struct simulation *sim)
{
	FILE* out;
	char temp[500];
	join_path(2,temp,get_output_path(sim),"log.dat");
	remove_file(sim,temp);
	//out=fopen(temp,"w");
	//fprintf(out,"gpvdm log file:\n");
	//fclose(out);

	join_path(2,temp,sim->root_simulation_path,"log_file_access.dat");
	remove_file(sim,temp);

	//out=fopen(temp,"w");
	//fprintf(out,"gpvdm file access log file:\n");
	//fclose(out);

}

void log_tell_use_where_file_access_log_is(struct simulation *sim)
{
	if (get_dump_status(sim,dump_file_access_log)==TRUE)
	{
		char temp[500];
		join_path(2,temp,sim->root_simulation_path,"log_file_access.dat");
		printf_log(sim,_("File access log written to %s\n"),temp);
	}

}
void log_write_file_access(struct simulation *sim,char * file,char mode)
{
	if (get_dump_status(sim,dump_file_access_log)==TRUE)
	{
		FILE* out;
		char temp[500];
		join_path(2,temp,sim->root_simulation_path,"log_file_access.dat");
		out=fopen(temp,"a");
		if (mode=='w')
		{
			fprintf(out,"write:%s\n",file);
		}else
		{
			fprintf(out,"read:%s\n",file);
		}

		fclose(out);
	}
}

void set_logging_level(struct simulation *sim,int value)
{
	sim->log_level=value;
}

void text_to_html(struct simulation *sim,char *out, char *in,int max)
{
	if (sim->html==TRUE)
	{

		int i=0;
		int out_pos=0;
		int len=0;
		len=strlen(in);
		for (i=0;i<len;i++)
		{
			if (in[i]=='\n')
			{
				out[out_pos]=0;
				strcat(out,"<br>");
				out_pos=strlen(out);
			}else
			{
				out[out_pos]=in[i];
				out_pos++;
				if (out_pos>=max)
				{
					out_pos--;
					break;
				}
			}

		}

		out[out_pos]=0;
	}else
	{
		out[0]=0;
		strncpy(out,in,max);
	}

}

void printf_log(struct simulation *sim, const char *format, ...)
{
	FILE* out;
	char data[STR_MAX];
	char data_html[STR_MAX];
	char temp[PATH_MAX];
	va_list args;
	va_start(args, format);
	vsnprintf(data,STR_MAX,format, args);

	if ((sim->log_level==log_level_screen)||(sim->log_level==log_level_screen_and_disk))
	{
		text_to_html(sim,data_html, data,1000);

		//wchar_t wide[1000];
		//int i=mbstowcs(wide, data_html, 1000);
		//wprintf(L"%S",wide);
			printf("%s",data_html);
			fflush(stdout);

	}

	if ((sim->log_level==log_level_disk)||(sim->log_level==log_level_screen_and_disk))
	{
		join_path(2,temp,get_output_path(sim),"log.dat");
		out=fopen(temp,"a");
		if (out==NULL)
		{
			wprintf(L"error: opening file %s\n",temp);
		}
		fprintf(out,"%s",data);
		fclose(out);
	}

	va_end(args);
}


void waveprint(struct simulation *sim,char *in,double wavelength)
{
	int r;
	int g;
	int b;

	if ((sim->log_level==log_level_screen)||(sim->log_level==log_level_screen_and_disk))
	{
		if (sim->html==TRUE)
		{
			wavelength_to_rgb(&r,&g,&b,wavelength*1e-9);
			printf_log(sim,"<font color=\"#%.2x%.2x%.2x\">",r,g,b);
		}else
		{
			if (wavelength<400.0)
			{
				textcolor(sim,fg_purple);
			}else
			if (wavelength<500.0)
			{
				textcolor(sim,fg_blue);
			}else
			if (wavelength<575.0)
			{
				textcolor(sim,fg_green);
			}else
			if (wavelength<600.0)
			{
				textcolor(sim,fg_yellow);
			}else
			{
				textcolor(sim,fg_red);

			}
		}
	}

	printf_log(sim,"%s",in);

	if ((sim->log_level==log_level_screen)||(sim->log_level==log_level_screen_and_disk))
	{
		if (sim->html==TRUE)
		{
			printf_log(sim,"</font>");
		}else
		{
			textcolor(sim,fg_reset);
		}

	}
}

void textcolor(struct simulation *sim, int color)
{
if (sim->html==TRUE)
{
	if (color==fg_purple)
	{
		printf_log(sim,"<font color=\"purple\">");
	}else
	if (color==fg_blue)
	{
		printf_log(sim,"<font color=\"blue\">");
	}else
	if (color==fg_green)
	{
		printf_log(sim,"<font color=\"green\">");
	}else
	if (color==fg_yellow)
	{
		printf_log(sim,"<font color=\"yellow\">");
	}else
	if (color==fg_red)
	{
		printf_log(sim,"<font color=\"red\">");
	}else
	if (color==fg_wight)
	{
		printf_log(sim,"<font color=\"#ffffff\">");
	}else
	if (color==fg_reset)
	{
		printf_log(sim,"</font>");
	}

}else
{
	printf("\033[%dm", color);
}

}

int log_search_error(char *path)
{
	int ret=-1;
    FILE * fp;
    char * line = NULL;
    ssize_t read;
	size_t len = 0;
    fp = fopen(path, "r");

    if (fp == NULL)
	{
       return ret;
	}

    while ((read = getline(&line, &len, fp)) != -1)
	{
		if (strcmp_begin(line,"error:")==0)
		{
			ret=0;
			break;
		}

    }

    fclose(fp);

    if (line)
	{
        free(line);
	}
    return ret;
}
